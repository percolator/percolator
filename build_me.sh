echo 'Creating out-of-source output folder...'
mkdir output
cd output
echo 'CMaking...'
cmake -DCMAKE_BUILD_TYPE=Debug ..
echo 'Building...'
sudo make install