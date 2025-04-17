#simple sh script that only works with having the original placement data
cd Z1907784975
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..
cd Z2732986066
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..
cd Z3343635604
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..
cd Z4324535763
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..
cd Z5185631889
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..
cd Z5185631911
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..
cd Z5348530222
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..
cd Z5348530626
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..
cd Z5348530683
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..
cd Z5348530741
cd placements
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
cd ..
mkdir lig_only
mv placements/lig* lig_only/
tar -czf lig_only.tar.gz lig_only/
cd ..

