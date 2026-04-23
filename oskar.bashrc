# ===== CUSTOM COMMANDS - OSKAR ===== #

# OSKAR Generic command
function oskar_bash() {
    local cflag=0
    local gflag=0
    local bflag=0
    local prog=""
    local ofile=""
    local prevd=$PWD
    local sifs=(${HOME}/.oskar/*.sif)
    local sfile="${sifs[0]}"

    if [[ ! -d ~/.oskar ]]; then
        mkdir ~/.oskar
        printf "${WARNING_TEXT}: OSKAR directory has just been created. This means that there is no existing OSKAR version on this device that can be recognised. Please download a binary or SIF file into ~/.oskar before continuing.\n"
    fi

    if [[ ! -f sfile ]]; then
        printf "${WARNING_TEXT}: There does not appear to be any SIF files in the ~/.oskar directory. Please make sure at least one sif file exists in the directory.\n"
    fi

    if [[ $1 == "--help" || $1 == "-h" ]]; then
        echo "=================================="
        echo "Run OSKAR with a custom generic bash command. By default it will run the first singularity image it"
        echo "finds in ~/.oskar. If there is no singularity image it will throw an error."
        echo "=================================="
        echo "Usage:"
        echo "oskar_bash [(-p|--prog)|(-s|--sif) <exec_file>] "
        echo "=================================="
        echo "Options:"
        echo "-g --global -s --sample Run OSKAR with sample settings"
        echo "-l --local              (Default) Run OSKAR in current directory with custom settings"
        echo "-i --intif              Run OSKAR's interferometer simulation"
        echo "-I --img                Run OSKAR's dirty imager simulation"
        echo "-b --beam               Run OSKAR's beam simulation"
        echo "-f --file               Settings file to use"
        echo "-c --clean              Clean directory of OSKAR logs"
        echo "-s --sif                If running singularity, what sif file to use"
        echo "-p --prog               If running a binary or application, the location of the binary or application"
        echo "=================================="

        return 0
    fi

    echo "Running custom OSKAR bash command ..."
    
    while [ $# -gt 0 ]; do
        case $1 in
            -g | --global | -s | --sample)
                gflag=1
            ;;
            -l | --local)
                gflag=0
            ;;
            -i | --intf)
                prog="oskar_sim_interferometer"
            ;;
            -b | --beam)
                prog="oskar_sim_beam_pattern"
            ;;
            -I | --image)
                prog="oskar_imager"
            ;;
            -f | --file)
                ofile=$2
                shift
            ;;
            -c | --clean)
                cflag=1
            ;;
            -s | --sif)
                sfile=$2
                bflag=0
                shift
            ;;
            -p | --prog)
                prog=$2
                bflag=1
                shift
            ;;
            \?)
                echo "'$1' is not a valid option. Use --help or -h to see what options are available."
            ;;
        esac
        shift
    done

    if [ $cflag -eq 1 ]; then
        if [ $gflag -eq 1 ]; then
            find ${HOME}/.oskar -name '*.log' -type f -delete
        else
            find . -name '*.log' -type f -delete
        fi
        return 0
    fi

    if [ $gflag -eq 1 ]; then
        ofile="$prog.ini"

        cd ${HOME}/.oskar
    fi

    if [ $bflag -eq 1 ]; then
        $prog $ofile
    else
        singularity exec --nv --bind $PWD --cleanenv --home $PWD $sfile $prog $ofile
    fi

    cd $prevd
}