cd build
cmake ..
make
stubgen -m komurasakiairbreathing
cd ..
mv build/out/komurasakiairbreathing.pyi komurasakiairbreathing.pyi