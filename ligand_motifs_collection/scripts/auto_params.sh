#Initially made by Matt Smith
#modified by Ari Ginsparg
source activate babel_env
module load python2
rawfile=$1
separate_path=$2
molf_to_params=$3

prot_code=$(echo $rawfile | sed 's/.pdb//g')
prot_chain=$(grep "^ATOM" $rawfile | head -n 1 |  awk '{ print $5 }')
lig_chain=$(grep "^HETATM" $rawfile | head -n 1 |  awk '{ print $5 }')
lig_name=$(grep "^HETATM" $rawfile | head -n 1 |  awk '{ print $4 }')


echo "File to process:"
echo $rawfile
echo "Ligand chain:"
echo $lig_chain
echo "Protein chain:"
echo $prot_chain
echo "Ligand name:"
echo $lig_name
echo "Protein code:"
echo $prot_code
echo $prot_code\_$prot_code

#./SeparatePDBsByChain.pl -pdbfile $rawfile -protchain $prot_chain > $prot_code\_$prot_code.pdb

#./SeparatePDBsByChain.pl -pdbfile $rawfile -ligchain $lig_chain | grep -v "HOH" > $rawfile.lig

${separate_path}SeparatePDBsByChain.pl -pdbfile $rawfile -protchain $prot_chain > $prot_code\_$prot_code.pdb

${separate_path}SeparatePDBsByChain.pl -pdbfile $rawfile -ligchain $lig_chain | grep -v "HOH" > $rawfile.lig


conda run babel -i pdb $rawfile.lig -o mol > $lig_name.mol



python ${molf_to_params}molfile_to_params.py --clobber $lig_name.mol -n $lig_name -p $prot_code

mv $lig_name.mol $prot_code

cat $prot_code\_0001.pdb >> $prot_code\_$prot_code.pdb
