source setup_stv.sh

PROCESSED_NTUPLE_DIR="/exp/uboone/data/users/barrow/CC2P_NewFiles/"
UNIV_OUTPUT_FILE=${PROCESSED_NTUPLE_DIR}"/Universes_CC2P.root"

MEASUREMENT_OUTPUT_FILE="./Output/"
UNF_MEAS_OUTPUT_FILE=${MEASUREMENT_OUTPUT_FILE}"/UnfoldedCrossSection_CC2P.root"

PELEE_NTUPLE_CONFIG="./Configs/files_to_process.txt"
FPM_CONFIG="./Configs/file_properties.txt"
BIN_CONFIG="./Configs/tutorial_bin_config.txt"
XSEC_CONFIG="./Configs/xsec_config.txt"
SLICE_CONFIG="./Configs/tutorial_slice_config.txt"
SYST_CONFIG="./Configs/systcalc.conf"

POSITIONAL_ARGS=()

DO_NTUPLE=0
DO_UNIVMAKE=0
DO_PLOT=0
DO_UNFOLD=0

while [[ $# -gt 0 ]]; do
  case $1 in
      -x|--xsec_config)
	  XSEC_CONFIG="$2"
	  shift # past argument
	  shift # past value
	  ;;
      -o|--output_file)
	  UNF_MEAS_OUTPUT_FILE="$2"
          shift # past argument
          shift # past value
          ;;
      -n|--ntuple)
	  DO_NTUPLE=1
	  shift # past argument
	  ;;
      -m|--univmake)
	  DO_UNIVMAKE=1
	  shift # past argument
	  ;;
      -p|--plot)
	  DO_PLOT=1
	  shift # past argument
	  ;;
      -u|--unfold)
	  DO_UNFOLD=1
	  shift # past argument
	  ;;
      -v|--validate)
	  DO_VALIDATE=1
	  shift # past argument
	  ;;
      -*|--*)
	  echo "Unknown option $1"
	  exit 1
	  ;;
      *)
	  POSITIONAL_ARGS+=("$1") # save positional arg
	  shift # past argument
	  ;;
  esac
done

if [[ $DO_NTUPLE -eq 1 ]]; then
    ./Scripts/ReprocessNTuples.sh ${PROCESSED_NTUPLE_DIR} ${PELEE_NTUPLE_CONFIG}
fi

if [[ $DO_UNIVMAKE -eq 1 ]]; then
    ./Scripts/UniverseMaker.sh ${FPM_CONFIG} ${BIN_CONFIG} ${UNIV_OUTPUT_FILE}
fi

if [[ $DO_PLOT -eq 1 ]]; then
    ./Scripts/PlotSlices.sh ${FPM_CONFIG} ${SYST_CONFIG} ${SLICE_CONFIG} ${UNIV_OUTPUT_FILE} ${MEASUREMENT_OUTPUT_FILE}
fi

if [[ $DO_UNFOLD -eq 1 ]]; then
    ./Scripts/Unfolder.sh ${XSEC_CONFIG} ${SLICE_CONFIG} ${UNF_MEAS_OUTPUT_FILE}
fi

if [[ $DO_VALIDATE -eq 1 ]]; then
    ./Scripts/Validate.sh ${XSEC_CONFIG} ${SLICE_CONFIG} ${UNF_MEAS_OUTPUT_FILE}
fi
