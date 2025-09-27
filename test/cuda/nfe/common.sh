check_AMBERHOME()
{
    if [ -z "${AMBERHOME}" ]; then
        echo " ** Error ** : AMBERHOME is not set"
        exit 1
    else
        if [ ! -d "${AMBERHOME}" ]; then
            echo " ** Error ** : \$AMBERHOME='${AMBERHOME}' is not a directory"
            exit 1
        fi
    fi
}

#
# sets SANDER variable
#

set_SANDER()
{
    if [ -z "${TESTsander}" ]; then
        check_AMBERHOME
        SANDER="${AMBERHOME}/bin/$1"
    else
        SANDER="${TESTsander}"
    fi

    if [ ! -x "${SANDER}" ]; then
        echo " ** Error ** : '${SANDER}' is not executable"
        exit 1
    fi
}

#
# sets SANDER_CMD (DO_PARALLEL SANDER) and NUMPROCS
#

set_SANDER_CMD()
{
    if [ "$1" -gt 1 ]; then
        if [ -z "${DO_PARALLEL}" ]; then
            echo " ** Error ** : empty DO_PARALLEL is not allowed in '`pwd`'"
            exit 1
        fi

        if [ ! -z "${MPICH_NP}" ]; then
            NUMPROCS="${MPICH_NP}"
        elif [ ! -z "${MP_PROCS}" ]; then
            NUMPROCS="${MP_PROCS}"
        else
            NUMPROCS=`echo ${DO_PARALLEL} | awk "BEGIN{i=0} {while (i < NF){if( \$ i == \"-np\" || \$ i == \"-n\"){i++;print \$ i; exit}; i++;}}" `
        fi

        if test -z "$NUMPROCS"; then
            echo " ** Error ** : could not determine NUMPROCS in '`pwd`'"
            exit 1
        fi

        NUMPROCS_is_okay="not so much"

        X2=`echo "$1 + $1" | bc`
        X4=`echo "$X2 + $X2" | bc`
        X8=`echo "$X4 + $X4" | bc`
        X16=`echo "$X8 + $X8" | bc`
        X32=`echo "$X16 + $X16" | bc`
        X64=`echo "$X32 + $X32" | bc`

        for x in $1 $X2 $X4 $X8 $X16 $X32 $X64; do
            if [ "${NUMPROCS}" -eq "$x" ]; then
                NUMPROCS_is_okay="yes it is"
                break
            fi
        done

        if [ "${NUMPROCS_is_okay}" != "yes it is" ]; then
            echo " This test case (`pwd`) requires a least $1 mpi threads."
            echo " The number of mpi threads must also be a multiple of $1"
            echo " and not more than ${X64}. Not running test, exiting ....."
            exit 1
        fi
    fi

    SANDER_CMD="${DO_PARALLEL} ${SANDER}"
}
#
# confirm DO_PARALLEL use fixed threads
#
set_SANDER_CMD_FIX()
{
    if [ "$1" -gt 1 ]; then
        if [ -z "${DO_PARALLEL}" ]; then
            echo " ** Error ** : empty DO_PARALLEL is not allowed in '`pwd`'"
            exit 1
        fi

        if [ ! -z "${MPICH_NP}" ]; then
            NUMPROCS="${MPICH_NP}"
        elif [ ! -z "${MP_PROCS}" ]; then
            NUMPROCS="${MP_PROCS}"
        else
            NUMPROCS=`echo ${DO_PARALLEL} | awk "BEGIN{i=0} {while (i < NF){if( \$ i == \"-np\" || \$ i == \"-n\"){i++;print \$ i; exit}; i++;}}" `
        fi

        if [ -z "$NUMPROCS" ]; then
            echo " ** Error ** : could not determine NUMPROCS in '`pwd`'"
            exit 1
        fi

        NUMPROCS_is_okay="not so much"

        X2=$1

        for x in $X2; do
            if [ "${NUMPROCS}" -eq "$x" ]; then
                NUMPROCS_is_okay="yes it is"
                break
            fi
        done

        if [ "${NUMPROCS_is_okay}" != "yes it is" ]; then
            echo " This test case (`pwd`) requires $X2 mpi threads."
            echo " Not running test, exiting ....."
            exit 1
        fi
    fi

    SANDER_CMD="${DO_PARALLEL} ${SANDER}"
}

set_NC_UTILS()
{
    check_AMBERHOME

    NCDUMP="${AMBERHOME}/bin/ncdump"
    NCGEN="${AMBERHOME}/bin/ncgen"

    if [ -x "${NCDUMP}" ]; then
        if test -x "${NCGEN}"; then
            return 0
        else
            unset NCDUMP NCGEN
            return 1
        fi
    else
        unset NCDUMP NCGEN
        return 1
    fi
}

do_ncdump()
{
    "${NCDUMP}" "$2" | sed -e "s|coeffs(row-major) ;|coeffs(row-major) ; coeffs:C_format=\"$1\" ;|g" | "${NCGEN}" -o ncdump.tmp ;
    "${NCDUMP}" ncdump.tmp > "$3" ; rm -f ncdump.tmp
}

save_junk_on_failure()
{
    FOLDER="junk.$$"
    DIFS=`find . -name \*.dif`
    if [ -n "${DIFS}" ]; then
        rm -rf "${FOLDER}"
        mkdir "${FOLDER}"
        for j in $*; do
            if [ -e $j ]; then
                cp -p $j "${FOLDER}"
            fi
        done
    fi
}

save_junk_on_failure()
{
    FOLDER="junk.$$"
    DIFS=`find . -name \*.dif`
    if [ -n "${DIFS}" ]; then
        rm -rf "${FOLDER}"
        mkdir "${FOLDER}"
        for j in $*; do
            if [ -e $j ]; then
                cp -p $j "${FOLDER}"
            fi
        done
    fi
}

CheckError() {
  if [ $1 -ne 0 ] ; then
    echo "${0}: Program error"
    exit 1
  fi
}

