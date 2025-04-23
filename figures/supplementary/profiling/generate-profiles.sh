mkdir profiles

# Done to ensure the scripts run from this directory
ln -s ../../Askja_VT-DLP_example/generate_results/inputs/ inputs

# Sweep over decimation factor, n_threads = 16, timestep = 300 s
memray run -o profiles/qm-detect.decimation-factor_1.n-threads_16.timestep_300.bin profile_detect.py --decimation_factor 1
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_16.timestep_300.bin profile_detect.py --decimation_factor 2
memray run -o profiles/qm-detect.decimation-factor_4.n-threads_16.timestep_300.bin profile_detect.py --decimation_factor 4

# Sweep over timestep, n_threads = 16, decimation-factor = 2
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_16.timestep_50.bin profile_detect.py --timestep 50
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_16.timestep_600.bin profile_detect.py --timestep 600
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_16.timestep_1200.bin profile_detect.py --timestep 1200
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_16.timestep_1800.bin profile_detect.py --timestep 1800

# Sweep over n_threads, decimation-factor = 2, timestep = 300 s
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_1.timestep_300.bin profile_detect.py --n_threads 1
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_2.timestep_300.bin profile_detect.py --n_threads 2
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_4.timestep_300.bin profile_detect.py --n_threads 4
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_8.timestep_300.bin profile_detect.py --n_threads 8
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_24.timestep_300.bin profile_detect.py --n_threads 24
memray run -o profiles/qm-detect.decimation-factor_2.n-threads_32.timestep_300.bin profile_detect.py --n_threads 32

# Now run each stage with the parameters used in the manuscript. Note, this may take some time!
memray run -o profiles/askja-lut.bin profile_lut.py
memray run -o profiles/askja-detect.bin profile_detect.py --n_threads 32 --timestep 300 --decimation_factor 2
memray run -o profiles/askja-trigger.bin profile_trigger.py
memray run -o profiles/askja-locate.bin profile_locate.py
