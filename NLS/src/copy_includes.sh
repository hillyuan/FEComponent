find . -name "*.h" -fprintf tmp.sh 'mkdir -p "../include/%h" \n cp "%h/%f" "../include/%h/"\n'
sed 's/\/.\//\//' tmp.sh > tmp_cp.sh
chmod a+x tmp_cp.sh
sh tmp_cp.sh
rm -rf tmp.sh tmp_cp.sh

