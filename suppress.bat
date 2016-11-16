cp -f sha/compressor_*.c .
make zipsrc
zip chaprosc sha/*
mv -f chaprosc.zip ../dist/shaprosc.zip
