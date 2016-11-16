unifdef -USHA -o compressor_prepare.c compressor_prepare.c
unifdef -USHA -o compressor_process.c compressor_process.c
make zipsrc
mv -f chaprosc.zip ../dist/chaprosc.zip
