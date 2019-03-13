# STV
Low energy starting track veto
STV.py is the main script. How to run:
python STV.py -i /data/sim/DeepCore/2018/pass2/gcd/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz /data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_000001_level12.i3.gz -o test --n
e 10 --sk 0 2>/dev/null

STV.py has an "STV_LE_Algorithm" icetray segment that you can import
STV_modules.py contains main functions
STV_utilities.py contains small supporting finction
STV_cuts.py has precuts
