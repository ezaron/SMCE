#Notes on the git-no-datafiles mirror on ebi:

cd SMCE
mkdir git-no-datafiles
cd git-no-datafiles
git clone https://github.com/ezaron/SMCE.git SMCE

cd SMCE
cp ../../fossil/edznc.jl .
cp ../../fossil/HAmod.jl .
cp ../../fossil/README.md .

mkdir HRET14 GDPdata RADSdata SWOTdata
cp ../../fossil/HRET14/README.md HRET14/
cp ../../fossil/HRET14/README.md GDPdata/
cp ../../fossil/HRET14/README.md RADSdata/
cp ../../fossil/HRET14/README.md SWOTdata/

