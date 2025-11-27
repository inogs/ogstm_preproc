#! /bin/bash

#SBATCH --job-name=Cl_CMEMS
#SBATCH -N1
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:20:00
#SBATCH --mem=300gb
#SBATCH --account=OGS_test2528
#SBATCH --partition=g100_meteo_prod
#SBATCH --qos=qos_meteo

module load autoload
module load intel/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary
module load mkl/oneapi-2021--binary
module load netcdf/4.7.4--oneapi--2021.2.0-ifort
module load netcdff/4.5.3--oneapi--2021.2.0-ifort
source /g100_work/OGS23_PRACE_IT/COPERNICUS/py_env_3.9.18_new/bin/activate

INPUTFILE="cmems_mod_med_bgc_my_4.2km-climatology_P1M-m_multi-vars_5.54W-36.29E_30.19N-45.98N_1.02-4152.90m_1999-01-01-1999-12-01.nc"
OUTDIR="_CLIM"


#if [ ! -f "$INPUTFILE" ]; then
#    python download_data.py
#else
#    echo "File $INPUTFILE already exists. Skipping download."
#fi

for VAR in o2_avg no3_avg chl_avg; do
	echo "processing var: $VAR"
	#python adjust_coordinates_dims.py -i "$INPUTFILE" -o "$OUTDIR" -v "$VAR"
done

OUTPUT_DIR=3D_LINKS
mkdir -p $OUTPUT_DIR

for YEAR in {2024..2025}; do
    # Ciclo sulle variabili CHL, O2, NO3
    for VAR in CHL O2 NO3; do
        # Assegna il nome corretto della variabile a vv
        if [[ "$VAR" == "CHL" ]]; then
            vv="P_l"
        elif [[ "$VAR" == "O2" ]]; then
            vv="O2o"
        else
            vv="N3n"
        fi
        
        INPUT_DIR="${VAR}_CLIM"  # Supponendo che le cartelle siano CHL_CLIM, O2_CLIM, NO3_CLIM
        
        # Controlla se la cartella esiste prima di procedere
        if [[ -d "$INPUT_DIR" ]]; then
            for file in "$INPUT_DIR"/ave.yyyy*15-00:00:00.${vv}.nc; do
                # Controlla se il file esiste per evitare errori
                [[ -e "$file" ]] || continue
                
                filename=$(basename "$file")

		MM=${filename:8:2}  # Estrai il mese dal nome del file
		#echo "Filename: $filename, Estratto MM: $MM"
                new_filename="ave.${YEAR}${MM}15-00:00:00.${vv}.nc"
		ln -s "$(realpath --relative-to="$OUTPUT_DIR" "$file")" "$OUTPUT_DIR/$new_filename"

            done
        fi
    done
done

#o2_std no3_std chl_std
