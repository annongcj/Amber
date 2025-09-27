"""
This is the program that contains the main updating facilities
"""
import os
from optparse import OptionParser, OptionGroup
import re
import shutil
import sys
from updateutils import downloader
from updateutils.downloader import HTTPError
from updateutils.exceptions import (IllegalOptions, UnknownArguments,
                                    NothingUpdated, SetupError)
from updateutils.preferences import Preferences
from updateutils.repos import LocalUnappliedPatchRepository as _LUPR
from updateutils.repos import LocalAppliedPatchRepository as _LAPR
from updateutils.repos import remote_patches
# Support Python 2.4 which doesn't have 'any'
if not 'any' in dir(__builtins__):
    from updateutils.utils import any
import warnings

if sys.version_info[0] < 3:
    input = raw_input

_PMEMD_VERSION = 24

# Until we know that our initial settings are OK (e.g., AMBERHOME is set
# correctly), do not allow any updates.

ALLOW_UPDATES = False


def _makes_changes(fcn):
    """
   Decorator that should be applied to any function in MainUpdater that actually
   applies or makes any other changes. This will bail out if ALLOW_UPDATES is
   False
   """

    def newfcn(self, *args, **kwargs):
        global ALLOW_UPDATES, WARNED
        if not ALLOW_UPDATES:
            raise SetupError(
                'Making changes has been disabled because PMEMDHOME '
                'is not set correctly or the updater has been moved from '
                'PMEMDHOME.')
        fcn(self, *args, **kwargs)

    return newfcn


class MainUpdater(object):
    """ This is the main class responsible for all updating activities """

    output = sys.stdout

    def __init__(self):
        """
      Sets up the remote and local repositories, loads the preferences, and does
      the necessary setup.
      """
        global _PMEMD_VERSION
        self.mustwarn = False
        self.preferences = Preferences('.updaterrc')
        self.verbose = None
        # Number of updates that are found via --check-updates
        self.num_updates = 0
        downloader.set_proxy_info(self.preferences['proxy'],
                                  self.preferences['proxy user'])
        if self.preferences['internet check'] is not None:
            downloader.DOMAIN_CHECK = self.preferences['internet check']

        # Set up the local bugfix repositories for Amber
        self.amber_unap = _LUPR(
            os.path.join('.patches',
                         'Amber%d_Unapplied_Patches' % _PMEMD_VERSION),
            'Amber %d' % _PMEMD_VERSION,
            prefix='update')
        self.amber_ap = _LAPR(
            os.path.join('.patches',
                         'Amber%d_Applied_Patches' % _PMEMD_VERSION),
            'Amber %d' % _PMEMD_VERSION,
            prefix='update')
        # Set the remote repositories. Start with the default names, then see if
        # the preferences defines a new location.
        amber = 'https://ambermd.org/bugfixes/pmemd/%d.0/' % _PMEMD_VERSION
        if self.preferences['amber updates'] is not None:
            amber = self.preferences['amber updates']
        self.amber_remote = remote_patches(
            amber, 'Amber %d' % _PMEMD_VERSION, 'update')
        # Set up a location for 3rd-party patches
        self.thirdparty_unap = _LUPR(
            os.path.join('.patches',
                         'ThirdParty%d_Unapplied_Patches' % _PMEMD_VERSION),
            '3rd Party Amber %d Patches' % _PMEMD_VERSION,
            prefix=None)
        self.thirdparty_ap = _LAPR(
            os.path.join('.patches',
                         'ThirdParty%d_Applied_Patches' % _PMEMD_VERSION),
            '3rd Party Amber %d Patches' % _PMEMD_VERSION,
            prefix=None)
        # Now set up all of the destination repositories (i.e., where remote
        # repos download to, where unapplied repos apply to, and where applied
        # repos reverse-apply to
        self.amber_remote.set_destination_repo(self.amber_unap)
        self.amber_unap.set_destination_repo(self.amber_ap)
        self.amber_ap.set_destination_repo(self.amber_unap)
        self.thirdparty_unap.set_destination_repo(self.thirdparty_ap)
        self.thirdparty_ap.set_destination_repo(self.thirdparty_unap)
        # Determine our patch level for both repositories
        self.amber_patchlevel = self.amber_ap.last()

    def version(self, dest=None, exit=None):
        """
      Determines the version of Amber we are running. Version number
      is reported as the major version -dot- bugfix number. If an exit code is
      given, sys.exit is called with that exit code. It will also print the
      string to a destination if dest is given
      """
        global _PMEMD_VERSION
        retstr = 'Version is reported as <version>.<patches applied>\n\n'
        retstr += '\n\t     Amber version %d.%s' % (
            _PMEMD_VERSION, str(self.amber_patchlevel).zfill(2))
        if not self.thirdparty_ap.empty():
            retstr += '\n\n\tThird-party modifications found. See full list of '
            retstr += '\n\tapplied patches for details.'
        if hasattr(dest, 'write'):
            dest.write(retstr + '\n\n')
        if isinstance(exit, int):
            sys.exit(exit)
        return retstr

    def set_timeout(self, value):
        """ Sets the default timeout """
        import socket
        socket.setdefaulttimeout(value)

    def check_updates(self):
        """ Looks for updates """
        global _PMEMD_VERSION
        # Now check for Amber updates
        self.output.write(
            '\nAvailable Amber %d patches:\n\n' % _PMEMD_VERSION +
            self.amber_remote.patch_details(self.amber_patchlevel + 1,
                                            self.verbose) + '\n')
        self.num_updates += self.amber_remote.last() - self.amber_ap.last()

    @_makes_changes
    def _update_tree(self, remote, localunap, localap, last):
        """
      Update the tree associated with the given remote, unapplied, and applied
      patch repositories
      """
        # Make sure there's something to do
        if last == 'all':
            last = remote.last()
        elif last <= localap.last():
            self.output.write(
                'Already at (or past) version %d. Cannot update\n' % (last))
            return
        elif last <= localunap.last():
            # Don't go online if we don't have to. Just apply from the unapplied
            # patch directory
            localunap.apply(first=localap.last() + 1, last=last)
            self.mustwarn = True
            return
        elif last > remote.last():
            return
        # Now print out a message if we aren't doing anything
        if last <= localap.last():
            self.output.write(
                'No new updates available for %s\n' % localap.reponame)
            return
        self.mustwarn = True
        remote.download_patches(
            first=localap.last() + 1, last=last, apply_=True)

    @_makes_changes
    def _update_amber(self, last='all'):
        """ Update Amber to the given patch level """
        self._update_tree(self.amber_remote, self.amber_unap,
                              self.amber_ap, last)

    @_makes_changes
    def _revert_tree(self, localap, target, desc):
        """ Reverse the desired tree to the target """
        # Make sure we have enough patches to rewind
        if localap.last() <= target:
            self.output.write('%s is at version %d. Not reverting to %d.\n' %
                              (desc, localap.last(), target))
            return
        if target < 0:
            self.output.write('Bad update number requested for %s (%d)\n' %
                              (desc, target))
            return
        self.mustwarn = True
        localap.apply(first=target + 1, last=localap.last())

    @_makes_changes
    def revert_amber(self, target):
        """ Reverses the Amber bugfix version """
        global _PMEMD_VERSION
        self._revert_tree(self.amber_ap, target,
                              'Amber %d' % _PMEMD_VERSION)

    def set_amber_remote(self, location):
        """ Sets a new amber remote location """
        global _PMEMD_VERSION
        self.preferences['amber updates'] = location
        self.amber_remote = remote_patches(
            location, 'Amber %d' % _PMEMD_VERSION, 'update')

    def finalize(self, normal=False):
        """ Finalizes the app -- saves preferences and such """
        self.preferences.save()
        self.thirdparty_unap.empty_repo()
        if normal and self.mustwarn:
            self.output.write(
                '\nNOTE: %s only updates the raw source code! You '
                'must recompile if you want \n      any changes to take effect!\n\n'
                % os.path.split(sys.argv[0])[1])

    def show_applied_patches(self):
        """ Show all applied patches """
        global _PMEMD_VERSION
        hdr = 'Amber %d Applied Patches:' % _PMEMD_VERSION
        self.output.write(hdr + '\n' + '-' * len(hdr) + '\n' +
                          self.amber_ap.patch_details(
                              verbose=self.verbose) + '\n\n')
        if not self.thirdparty_ap.empty():
            hdr = 'Third-party Patches:'
            self.output.write(hdr + '\n' + '-' * len(hdr) + '\n' +
                              self.thirdparty_ap.patch_details(
                                  verbose=self.verbose) + '\n\n')

    def show_unapplied_patches(self):
        """ Shows all unapplied patches that have been downloaded """
        global _PMEMD_VERSION
        self.output.write(
            '\nAmber %d Unapplied Patches:\n\n' % _PMEMD_VERSION +
            self.amber_unap.patch_details(verbose=self.verbose) + '\n')
        if not self.thirdparty_unap.empty():
            self.output.write(
                'Third-party Unapplied Patches:\n',
                self.thirdparty_unap.patch_details(verbose=self.verbose))

    @_makes_changes
    def download_patches(self):
        """
      Downloads all available patches that have not yet been downloaded or
      applied
      """
        self.output.write(
            'Looking for and downloading updates... please wait.\n')
        self.amber_remote.download_patches(
            self.amber_ap.last() + 1, apply_=False)

    @_makes_changes
    def update(self):
        """ Updates to the current version with all patches applied """
        self.output.write('Preparing to apply updates... please wait.\n')
        self._update_amber()

    @_makes_changes
    def update_to(self, version):
        """ Updates to a given version based on a passed VERSION string """
        self.output.write(
            'Preparing to apply specific updates... please wait.\n')
        update_applied = False
        myre = re.compile(r'amber[/\.](\d+)', re.IGNORECASE)
        updateto = myre.findall(version)
        if updateto:
            self._update_amber(int(updateto[0]))
            update_applied = True
        if not update_applied:
            warnings.warn('--update-to: no matching VERSION found',
                          NothingUpdated)

    @_makes_changes
    def revert_to(self, version):
        """ Reverts to a given version based on a passed VERSION string """
        self.output.write(
            'Preparing to undo specific updates... please wait.\n')
        reverted = False
        myre = re.compile(r'amber[/\.](\d+)', re.IGNORECASE)
        revertto = myre.findall(version)
        if revertto:
            self.revert_amber(int(revertto[0]))
            reverted = True
        if not reverted:
            warnings.warn('--revert-to: no matching VERSION found',
                          NothingUpdated)

    def get_user_input(self, message):
        """ Gets input from the user """
        return input(message.strip() + ' ')

    @_makes_changes
    def thirdparty_apply(self, patch):
        """ Applies a third-party patch """
        # If the patch is online, download it. Otherwise, copy it
        self.output.write('Preparing to apply a third-party update/patch... '
                          'please wait.\n')
        if patch.startswith('http://'):
            try:
                dlmger = downloader.DownloadManager(patch,
                                                    os.path.split(patch)[1])
                dlmger.download_to(self.thirdparty_unap.location +
                                   os.path.split(patch)[1])
            except HTTPError:
                self.output.write(
                    '%s could not be downloaded... Quitting...\n' % patch)
                sys.exit(1)
        else:
            try:
                shutil.copy(
                    patch,
                    self.thirdparty_unap.location + os.path.split(patch)[1])
            except IOError:
                self.output.write(
                    '%s could not be found/copied... Quitting...\n' % patch)
                sys.exit(1)
        # Reload the third-party repo and apply them
        self.thirdparty_unap.fill_list()
        self.thirdparty_unap.apply(self.thirdparty_unap.location +
                                   os.path.split(patch)[1])
        self.mustwarn = True

    @_makes_changes
    def thirdparty_reverse(self, patch):
        """ Reverses a third-party patch that has been applied """
        try:
            idx = self.thirdparty_ap.patch_list.index(patch)
            self.thirdparty_ap.apply(first=idx, last=idx)
            self.mustwarn = True
        except ValueError:
            self.output.write(
                ('%s could not be found in the applied patches... '
                 'Quitting...\n') % patch)

    @_makes_changes
    def remove_unapplied(self):
        """ Removes all unapplied patches """
        self.amber_unap.empty_repo()
        self.thirdparty_unap.empty_repo()

    def mainloop(self):
        """ Executes the commands in the patcher """
        global _PMEMD_VERSION
        self.define_clparser()
        # This call is basically just used to direct --help to the basic options,
        # only. The real parsing is done by the full parser
        opt, arg = self.full_parser.parse_args()
        # Check that all options are OK
        try:
            self._check_options(opt, arg)
        except IllegalOptions:
            err = sys.exc_info()[1]
            self.output.write('IllegalOptions: %s\n' % err)
            sys.exit(1)
        if len(sys.argv) == 1:
            self.basic_parser.print_help()
            sys.exit(1)
        # First go through and update the preferences
        if opt.proxy is not None:
            self.preferences['proxy'] = opt.proxy
        if opt.proxy_user is not None:
            self.preferences['proxy user'] = opt.proxy_user
        if opt.delete_proxy:
            self.preferences['proxy'] = None
            self.preferences['proxy user'] = None
        if opt.amber_updates is not None:
            self.set_amber_remote(opt.amber_updates)
        if opt.reset_remotes:
            self.set_amber_remote(
                'https://ambermd.org/bugfixes/pmemd/%d.0/' % _PMEMD_VERSION)
            self.preferences['amber updates'] = None
        # Now set the verbose setting and the timeout value
        self.verbose = opt.verbose
        self.set_timeout(opt.timeout)
        # Now do any querying options
        if opt.check:
            self.output.write(
                'Checking for available patches online. This may '
                'take a few seconds...\n')
            self.check_updates()
        if opt.remove_unapplied:
            self.remove_unapplied()
        if opt.show_applied:
            self.show_applied_patches()
        if opt.show_unapplied:
            self.show_unapplied_patches()
        if opt.download:
            self.download_patches()
        if opt.update:
            self.update()
        if opt.update_to:
            self.update_to(opt.update_to)
        if opt.revert_to:
            self.revert_to(opt.revert_to)
        if opt.apply:
            self.thirdparty_apply(opt.apply)
        if opt.reverse:
            self.thirdparty_reverse(opt.reverse)

    def _check_options(self, opt, arg):
        """ Checks for illegal option combinations """
        if arg:
            warnings.warn('Unknown arguments are being ignored: %s' % arg,
                          UnknownArguments)
        if opt.verbose is not None and (opt.verbose < 0 or opt.verbose > 4):
            raise IllegalOptions('Verbose must be 0, 1, 2, 3, or 4')
        if opt.timeout <= 0:
            raise IllegalOptions('timeout must be greater than zero')
        if opt.update and any((opt.update_to, opt.revert_to)):
            raise IllegalOptions(
                'All Version Control update options are mutually '
                'exclusive')
        if opt.update_to and any((opt.update, opt.revert_to)):
            raise IllegalOptions(
                'All Version Control update options are mutually '
                'exclusive')
        if opt.revert_to and any((opt.update, opt.update_to)):
            raise IllegalOptions(
                'All Version Control update options are mutually '
                'exclusive')
        if opt.apply and opt.reverse:
            raise IllegalOptions(
                'Third-party patches must be applied/reversed '
                'independently.')
        if opt.amber_updates and opt.reset_remotes:
            raise IllegalOptions('Cannot set and delete update repositories '
                                 'simultaneously')
        if any((opt.proxy, opt.proxy_user)) and opt.delete_proxy:
            raise IllegalOptions('Cannot set up and delete proxy settings '
                                 'simultaneously')

    def define_clparser(self):
        """ Defines the command-line parser """
        global _PMEMD_VERSION
        from updateutils.utils import print_updater_version
        # Set up a basic parser to print easy-going help
        try:
            parser = OptionParser(epilog='This program handles bug fixes and '
                                  'updates released for Amber.')
        except TypeError:
            parser = OptionParser()
        parser.add_option(
            '--full-help',
            action='callback',
            help='Print the help '
            'for all available options, including advanced options.',
            callback=lambda *args, **kwargs: 0)
        parser.add_option(
            '-v',
            '--version',
            action='callback',
            help='Display Amber version.',
            callback=lambda *args, **kwargs: self.version(self.output, 0))
        parser.add_option('-V', action='callback',
              help='Print version of updater script.',
              callback=lambda arg, opt, value, parser:
                 print_updater_version(self, self.output, 0))
        group = OptionGroup(parser, 'Update Options', 'These are the basic '
                            'options for updating Amber')
        group.add_option(
            '--check-updates',
            dest='check',
            action='store_true',
            default=False,
            help='Check online for any available updates.')
        group.add_option(
            '--show-applied-patches',
            dest='show_applied',
            action='store_true',
            help='Show details of all patches that have '
            'been applied so far.')
        group.add_option(
            '--update',
            dest='update',
            action='store_true',
            default=False,
            help='Update to the latest version of '
            'Amber (apply all available updates).')
        parser.add_option_group(group)
        group = OptionGroup(parser, 'Connection Settings',
                            'These options allow '
                            'you to customize how the internet is accessed.')
        group.add_option(
            '--proxy',
            dest='proxy',
            default=None,
            metavar='<PROXY ADDRESS>',
            help='All connections to the internet '
            'goes through the specified proxy.')
        group.add_option(
            '--proxy-user',
            dest='proxy_user',
            default=None,
            metavar='<USERNAME>',
            help='If your proxy is authenticated, this is '
            'the username that will be used in authentication. You will be '
            'prompted for your password the first time internet access is '
            'needed.')
        group.add_option(
            '--delete-proxy',
            dest='delete_proxy',
            default=False,
            action='store_true',
            help='Delete stored proxy information.')
        parser.add_option_group(group)
        self.basic_parser = parser

        try:
            parser = OptionParser(
                epilog='This program handles bug fixes and '
                'updates released for Amber.',
                add_help_option=False)
        except TypeError:
            parser = OptionParser(add_help_option=False)
        parser.add_option(
            '-h',
            '--help',
            action='callback',
            help='show basic '
            'help message and exit.',
            callback=
            lambda *args, **kwargs: self.basic_parser.parse_args(['-h']))
        parser.add_option(
            '--full-help',
            action='callback',
            help='show this '
            'help message and exit.',
            callback=
            lambda var, opt, value, parser: sys.exit(parser.print_help()))
        parser.add_option(
            '-v',
            '--version',
            action='callback',
            help='Display Amber version',
            callback=lambda *args, **kwargs: self.version(self.output, 0))
        parser.add_option('-V', action='callback',
              help='Print version of updater script.',
              callback=lambda arg, opt, value, parser:
                 print_updater_version(self, self.output, 0))
        parser.add_option(
            '--verbose',
            dest='verbose',
            type='int',
            default=None,
            metavar='<0|1|2|3|4>',
            help='Controls the amount of information '
            'printed about all patches. Larger values print more information.')
        parser.add_option(
            '--timeout',
            dest='timeout',
            type='float',
            metavar='<SECONDS>',
            default=10.0,
            help='How long to wait for a '
            'response from the remote repository. Default is %default seconds.')
        group = OptionGroup(
            parser, 'Querying Options', 'These options are used '
            'to obtain information about patches that have been applied, '
            'are available, and that have been downloaded but not yet '
            'applied.')
        group.add_option(
            '--check-updates',
            dest='check',
            action='store_true',
            default=False,
            help='Check online for any available updates.')
        group.add_option(
            '--show-applied-patches',
            dest='show_applied',
            action='store_true',
            help='Show details of all patches that have '
            'been applied so far.')
        group.add_option(
            '--show-unapplied-patches',
            dest='show_unapplied',
            action='store_true',
            help='Show details of all patches that have '
            'been downloaded but not applied yet.')
        parser.add_option_group(group)
        group = OptionGroup(
            parser, 'Version Control', 'These options are used '
            'to control which bug fixes and updates are applied and which '
            'may be reversed. VERSION strings are comma-separated values '
            'Amber/<bugfix> (case-insensitive).')
        group.add_option(
            '--update',
            dest='update',
            action='store_true',
            default=False,
            help='Update to the latest version of '
            'Amber (apply all available updates).')
        group.add_option(
            '--update-to',
            dest='update_to',
            metavar='<VERSION>',
            default=None,
            help='Apply as many updates as necessary to get to '
            'the desired "dot" version(s) of Amber. Will NOT '
            'reverse any patches.'
            '   Example: --update-to=Amber.3')
        group.add_option(
            '--revert-to',
            dest='revert_to',
            metavar='<VERSION>',
            default=None,
            help='Reverse as many updates as necessary to get to '
            'the desired "dot" version(s) of Amber. '
            'Will NOT apply any new patches.'
            '   Example: --revert-to=Amber.3')
        group.add_option(
            '--download-patches',
            dest='download',
            default=False,
            action='store_true',
            help='Download available patches, but do not '
            'apply them.')
        group.add_option(
            '--remove-unapplied',
            dest='remove_unapplied',
            default=False,
            action='store_true',
            help='Removes all patches that '
            'have previously been downloaded, but not yet applied.')
        parser.add_option_group(group)
        group = OptionGroup(
            parser, 'Third-party Patches', 'These options '
            'control applying third-party patches to customize '
            'Amber from sources besides Amber.')
        group.add_option(
            '--apply-patch',
            dest='apply',
            default=None,
            metavar='<PATCH>',
            help='Apply the provided third-party patch. If '
            'the patch starts with http://, it will be downloaded from the '
            'web and applied. If it is local, it will be applied as a '
            'third-party patch.')
        group.add_option(
            '--reverse-patch',
            dest='reverse',
            default=None,
            metavar='<PATCH>',
            help='Reverse the given patch name. This patch '
            'name must have been a previously applied, third-party patch. If '
            'present, it will be reversed.')
        parser.add_option_group(group)
        group = OptionGroup(
            parser, 'Connection Settings', 'These options allow '
            'you to customize how online resources will be accessed and '
            'patches downloaded (if you wish to set up any mirrors).')
        group.add_option(
            '--amber-updates',
            dest='amber_updates',
            default=None,
            metavar='<URL | FOLDER>',
            help='New location to search for updates '
            'to Amber %d' % _PMEMD_VERSION)
        group.add_option(
            '--proxy',
            dest='proxy',
            default=None,
            metavar='<PROXY ADDRESS>',
            help='All connections to the internet '
            'goes through the specified proxy.')
        group.add_option(
            '--proxy-user',
            dest='proxy_user',
            default=None,
            metavar='<USERNAME>',
            help='If your proxy is authenticated, this is '
            'the username that will be used in authentication. You will be '
            'prompted for your password the first time internet access is '
            'needed.')
        group.add_option(
            '--internet-check',
            dest='internet_check',
            default=None,
            metavar='<URL>',
            help='Website to check for access prior to '
            'downloading patches. Defaults to https://ambermd.org')
        group.add_option(
            '--delete-proxy',
            dest='delete_proxy',
            default=False,
            action='store_true',
            help='Delete stored proxy information.')
        group.add_option(
            '--reset-remotes',
            dest='reset_remotes',
            default=False,
            action='store_true',
            help='Restore online update sources to their '
            'original defaults.')
        group.add_option(
            '--reset-check',
            dest='reset_check',
            default=False,
            action='store_true',
            help='Reset website to test for internet '
            'access to https://ambermd.org')
        parser.add_option_group(group)
        self.full_parser = parser
