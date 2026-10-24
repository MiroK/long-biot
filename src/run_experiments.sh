nrefs=5

# #Show failure is a combo
# pbcs="PPFF"
# python3 biot2.py -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 80 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs

# # It does not happen with other bcs
# pbcs="PPPP"
# python3 biot2.py -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
# python3 biot2.py -L 80 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs

# # Or when c is large enough
# pbcs="PPFF"
# python3 biot2.py -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
# python3 biot2.py -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
# python3 biot2.py -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
# python3 biot2.py -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
# python3 biot2.py -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs
# python3 biot2.py -L 80 -precond standard -ubcs TTDD -pbcs $pbcs -c 1.0 -nrefs $nrefs

# #Fail of other formulations
# for script in biot3.py biot4.py
# do
#     pbcs="PPFF"
#     python3 $script -L 1 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
#     python3 $script -L 5 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
#     python3 $script -L 10 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
#     python3 $script -L 20 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
#     python3 $script -L 40 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs
#     python3 $script -L 80 -precond standard -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs    
# done
    
# Fix by espen [exact]
for script in biot2.py 
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


# for script in biot3.py 
# do
#     pbcs="PPFF"
#     python3 $script -L 1 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 5 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 10 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 20 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 40 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 80 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg    
# done

# for script in biot4.py 
# do
#     pbcs="PPFF"
#     python3 $script -L 1 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 5 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 10 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 20 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 40 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg
#     python3 $script -L 80 -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ pyamg    
# done


# <<<<<<< Updated upstream
# nrefs=5
# pbcs="PPFF"
# for L in 1 20 40 80
# do
#     for script in biot3.py biot4.py
#     do
# 	# Loop lmbda
# 	for lmbda in 1E10 1E8 1E0 1E2 1E4 1E6
# 	do
# 	    for K in 1E-10 1E0 1E-2 1E-4 1E-6 1E-8
# 	    do
# 		python3 $script -L $L -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ lu -lmbda $lmbda -K $K
# =======
# nrefs=3
# pbcs="PPFF"
# for L in 2
# do
#     for script in biot3.py biot4.py
#     do
# 	for K in 1E0 1E-2 1E-4 1E-6 1E-8 1E-10
# 	do
# 	    for lmbda in 1E0 1E2 1E4 1E6 1E8 1E10
# 	    do
# 		python3 $script -L $L -precond espen -ubcs TTDD -pbcs $pbcs -c 0.0 -nrefs $nrefs -inverseQ lu -K $K -lmbda $lmbda
# >>>>>>> Stashed changes
# 	    done
# 	done
#     done
# done
