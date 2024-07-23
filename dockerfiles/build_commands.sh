#Using docker version 24.0.5 and singularity version 4.1.0
#$1 is the name of the container that you want, which will pass down all steps

#build container:
echo Making container with name: $1
sudo docker build - < Dockerfile -t $1 --progress=plain

#get image id:
image_id=$(sudo docker images | grep $1 | awk '{print $3}')
echo $1 image id: $image_id

#write image to tar that singularity can use:
echo Saving image to tar
sudo docker save $image_id -o $1.tar

#build sif container from image
echo Converting tar to Singularity SIF: $1.sif
sudo singularity build $1.sif docker-archive://$1.tar
