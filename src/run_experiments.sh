nrefs=5
#Show failure is a combo
pbcs="PPFF"
python3 biot2.py -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 80 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs

# It does not happen with other bcs
pbcs="PPPP"
python3 biot2.py -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
python3 biot2.py -L 80 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs

# Or when c is large enough
pbcs="PPFF"
python3 biot2.py -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
python3 biot2.py -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
python3 biot2.py -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
python3 biot2.py -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
python3 biot2.py -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
python3 biot2.py -L 80 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs

#Fail of other formulations
for script in biot3.py biot4.py
do
    pbcs="PPFF"
    python3 $script -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 80 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs    
done
    
# Fix by espen [exact]
for script in biot2.py biot3.py biot4.py
do
    pbcs="PPFF"
    python3 $script -L 1 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 5 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 10 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 20 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 40 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
    python3 $script -L 80 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs    
done

# Approximate
for script in biot2.py 
do
    pbcs="PPFF"
    python3 $script -L 1 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ amg
    python3 $script -L 5 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ amg
    python3 $script -L 10 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ amg
    python3 $script -L 20 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ amg
    python3 $script -L 40 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ amg
    python3 $script -L 80 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ amg    
done


for script in biot3.py 
do
    pbcs="PPFF"
    python3 $script -L 1 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 5 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 10 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 20 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 40 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 80 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg    
done

for script in biot4.py 
do
    pbcs="PPFF"
    python3 $script -L 1 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 5 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 10 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 20 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 40 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
    python3 $script -L 80 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg    
done
