#!/bin/bash
set -e # exit on errors
set -x # print commands and their arguments as they are executed

if [ -e VARS2 ]; then
  source VARS2
fi

find_normaliz_binary() {
    if [ -x "${PWD}/source/normaliz" ]; then
        echo "${PWD}/source/normaliz"
    elif [ -x "${PWD}/../source/normaliz" ]; then
        echo "${PWD}/../source/normaliz"
    else
        echo ""
    fi
}

print_failed_diffs() {
    local test_dir="$1"

    if [ ! -d "${test_dir}" ]; then
        return
    fi

    echo "==== Non-empty test diffs ===="
    find "${test_dir}" -name '*.diff' -type f -size +0 -print0 |
        while IFS= read -r -d '' diff_file; do
            echo "---- ${diff_file} ----"
            sed -n '1,160p' "${diff_file}"
        done
}

rerun_sigill_probe() {
    local normaliz="$1"
    local test_dir="$2"
    local input="$3"

    if [ ! -x "${normaliz}" ] || [ ! -d "${test_dir}" ]; then
        return
    fi

    if [ ! -f "${test_dir}/${input}.in" ]; then
        return
    fi

    echo "==== SIGILL probe: ${input} ===="
    (
        set +e
        cd "${test_dir}" || exit 0
        if /usr/bin/time -v true > /dev/null 2>&1; then
            /usr/bin/time -v "${normaliz}" -c --inv "${input}" > /dev/null
        else
            time "${normaliz}" -c --inv "${input}" > /dev/null
        fi
        local status=$?
        echo "probe exit status: ${status}"

        if ! command -v gdb > /dev/null 2>&1; then
            echo "gdb not available; trying to install it"
            if command -v apt-get > /dev/null 2>&1; then
                if command -v sudo > /dev/null 2>&1; then
                    sudo apt-get update
                    sudo apt-get install -y gdb
                else
                    apt-get update
                    apt-get install -y gdb
                fi
            fi
        fi

        if command -v gdb > /dev/null 2>&1; then
            gdb -batch \
                -ex 'run > /dev/null' \
                -ex 'x/i $pc' \
                -ex bt \
                --args "${normaliz}" -c --inv "${input}"
        else
            echo "gdb not available"
        fi
    )
}

diagnose_check_failure() {
    local status="$1"
    local normaliz
    local test_dir="test/run_tests"

    set +e
    echo "==== make check failed with status ${status} ===="
    echo "==== Host CPU ===="
    if command -v lscpu > /dev/null 2>&1; then
        lscpu
    else
        uname -a
        sysctl -a 2>/dev/null | grep -E 'machdep.cpu|hw.optional' || true
    fi

    normaliz="$(find_normaliz_binary)"
    if [ -n "${normaliz}" ]; then
        echo "==== normaliz binary ===="
        file "${normaliz}"
        if [[ $OSTYPE == darwin* ]]; then
            otool -L "${normaliz}"
        else
            ldd "${normaliz}"
        fi
    else
        echo "normaliz binary not found"
    fi

    print_failed_diffs "${test_dir}"

    rerun_sigill_probe "${normaliz}" "${test_dir}" "test-/rational_inhom_H"
    rerun_sigill_probe "${normaliz}" "${test_dir}" "test-/rational_inhom_full"
    rerun_sigill_probe "${normaliz}" "${test_dir}" "test-/rational_inhom_hedron"

    exit "${status}"
}

case $BUILDSYSTEM in

    *static*)
        ## export NORMPARA=-x=1 ## if paralleization should not work

        make install
        if [[ $OSTYPE == darwin* ]]; then
            echo "WWWWWWWWWWWWWWWWWWWWWW"
            diff --version
            echo "VVVVVVVVVVVVVVVVVVVVVV"
            otool -L ${PREFIX}/bin/*
            echo "UUUUUUUUUUUUUUUUUUUUUU"
        else
            ldd ${PREFIX}/bin/*
        fi

        make check

        cd local/bin
        echo "CONTENTS OF LOCAL/BIN"
        ls .
        if [ "x$INTEL" != x ]; then
            zip MacOSbinary-Intel.zip *
        else
            zip MacOSbinary.zip *
        fi

        ;;

    *makedistcheck*)
        make -j2 distcheck
        ;;

    *)
        # make -j2 -k
        make -j2 -k check || diagnose_check_failure $?
        make install
        if [[ $OSTYPE == darwin* ]]; then
            echo "VVVVVVVVVVVVVVVVVVVVVV"
            diff --version
            echo "VVVVVVVVVVVVVVVVVVVVVV"
            otool -L ${PREFIX}/bin/*
        else
            ldd ${PREFIX}/bin/*
        fi
        make installcheck
        ;;
esac

export -p | sed 's/declare -x/export/g' > VARS3
