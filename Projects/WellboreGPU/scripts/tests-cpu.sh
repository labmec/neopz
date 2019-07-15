#!/bin/bash
rm -rf CPU

mkdir CPU

echo Configuring compilation options with CMake
echo Test type:		CPU
echo Polynomial order:	1	
cmake -DCMAKE_BUILD_TYPE=Release -DUSING_CUDA=off -DUSING_TBB=on -DUSING_SPARSE=on -DCOMPUTE_K_HYBRID=off -DO_LINEAR=on -DO_QUADRATIC=off -DO_CUBIC=off . > any.txt

echo Compiling...
make -j32 > any.txt

for j in 3 15 64 260 1044
do
    echo "\nMesh:		$j k"
    mkdir CPU/order1-mesh$j
    for i in 1 2 3 4 5
    do
        echo $i of 5
        ./IntPointsFEM $j > CPU/order1-mesh$j/out-$i.txt
    done
done

echo ""
echo Configuring compilation options with CMake
echo Test type:         CPU
echo Polynomial order:  2       
cmake -DCMAKE_BUILD_TYPE=Release -DUSING_CUDA=off -DUSING_TBB=on -DUSING_SPARSE=on -DCOMPUTE_K_HYBRID=off -DO_LINEAR=off -DO_QUADRATIC=on -DO_CUBIC=off . > any.txt

echo Compiling...
make -j32 > any.txt

for j in 3 15 64 260 1044
do
    echo "\nMesh:               $j k"
    mkdir CPU/order2-mesh$j
    for i in 1 2 3 4 5
    do
        echo $i of 5
        ./IntPointsFEM $j > CPU/order2-mesh$j/out-$i.txt
    done
done


echo ""
echo Configuring compilation options with CMake
echo Test type:         CPU
echo Polynomial order:  3       
cmake -DCMAKE_BUILD_TYPE=Release -DUSING_CUDA=off -DUSING_TBB=on -DUSING_SPARSE=on -DCOMPUTE_K_HYBRID=off -DO_LINEAR=off -DO_QUADRATIC=off -DO_CUBIC=on . > any.txt

echo Compiling...
make -j32 > any.txt

for j in 3 15 64 260
do
    echo "\nMesh:               $j k"
    mkdir CPU/order3-mesh$j
    for i in 1 2 3 4 5
    do
        echo $i of 5
        ./IntPointsFEM $j > CPU/order3-mesh$j/out-$i.txt
    done
done



rm any.txt



