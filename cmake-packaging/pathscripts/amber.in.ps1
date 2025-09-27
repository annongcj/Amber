# Run this script to add the variables necessary to use Amber to your shell.
# This script must be located in the Amber root folder!

# To run this script, you must change your Powershell execution policy to unrestricted, with the command "Set-ExecutionPolicy Unrestricted".

# Amber was configured on @CONFIGURATION_TIMESTAMP@
	
$env:AMBERHOME= Split-Path -Path $MyInvocation.MyCommand.Path
$env:PATH= ";" + $env:AMBERHOME + "\bin;" + "@EXTRA_PATH_FOLDERS_PSH@" + $env:PATH

# add Perl libs to path (only enabled if Amber was compiled with Perl support)
If(@AMBERSH_PERL_SUPPORT@ -eq 1) {
	If([string]::IsNullOrEmpty($env:PERL5LIB)) {
		$env:PERL5LIB= $env:AMBERHOME + "\@PERL_MODULE_DIR_WIN@"
	}
	else {
		$env:PERL5LIB= ";" + $env:AMBERHOME + "\@PERL_MODULE_DIR_WIN@" + $env:PERL5LIB
	}
}

# add Python libs to path (only enabled if Amber was compiled with Python support)
If(@AMBERSH_PYTHON_SUPPORT@ -eq 1) {
	If([string]::IsNullOrEmpty($env:PYTHONPATH)) {
		$env:PYTHONPATH= $env:AMBERHOME + "@RELATIVE_PYTHONPATH_WIN@"
	}
	else {
		$env:PYTHONPATH= ";" + $env:AMBERHOME + "@RELATIVE_PYTHONPATH_WIN@" + $env:PYTHONPATH
	}
}

# Tell QUICK where to find its data files
If(@AMBERSH_QUICK_SUPPORT@ -eq 1) {
	$env:PERL5LIB= $env:AMBERHOME + "\AmberTools\src\quick\basis"
}
