# Show failure is a combo
pbcs="PPFF"
python3 biot2.py -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6
python3 biot2.py -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6
python3 biot2.py -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5
python3 biot2.py -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5
python3 biot2.py -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 4

# It does not happen with other bcs
pbcs="PPPP"
python3 biot2.py -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6
python3 biot2.py -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6
python3 biot2.py -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5
python3 biot2.py -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5
python3 biot2.py -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 4

# Or when c is large enough
pbcs="PPFF"
python3 biot2.py -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 1 -nrefs 6
python3 biot2.py -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 1 -nrefs 6
python3 biot2.py -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 1 -nrefs 5
python3 biot2.py -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 1 -nrefs 5
python3 biot2.py -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 1 -nrefs 4

# Fail of other formulations
for script in biot3.py biot4.py
do
    pbcs="PPFF"
    python3 $script -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6
    python3 $script -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6
    python3 $script -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5
    python3 $script -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5
    python3 $script -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 4
done
    
# Fix by espen [exact]
for script in biot2.py biot3.py biot4.py
do
    pbcs="PPFF"
    python3 $script -L 1 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6
    python3 $script -L 5 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6
    python3 $script -L 10 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5
    python3 $script -L 20 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5
    python3 $script -L 40 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 4
done

# Approximate
for script in biot2.py biot3.py biot4.py
do
    pbcs="PPFF"
    python3 $script -L 1 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6 -inverseQ amg
    python3 $script -L 5 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 6 -inverseQ amg
    python3 $script -L 10 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5 -inverseQ amg
    python3 $script -L 20 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 5 -inverseQ amg
    python3 $script -L 40 -precond espen -ubcs TTDD -pbcs $pbcs -c 0 -nrefs 4 -inverseQ amg
done
