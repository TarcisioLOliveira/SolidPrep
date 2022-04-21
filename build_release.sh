mkdir -p release
cd release
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release .. -DCMAKE_INSTALL_PREFIX=. && make install
