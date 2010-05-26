echo 'Creating out-of-source output folder...'
mkdir output
cd output
echo 'CMaking...'
cmake -DCMAKE_BUILD_TYPE=Release -DGOOGLE_TEST=FALSE ..
echo 'Building...'
make
echo 'Installing...'
sudo make install