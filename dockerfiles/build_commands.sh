#Using docker version 24.0.5 and singularity version 4.1.0

#build shapedb container:
sudo docker build - < Dockerfile -t shapedb_container:1.2.1 --progress=plain

#get image id:

#write image to tar that singularity can use:
sudo docker save $IMAGE_ID -o shapedb_container.tar

#build sif container from image
sudo singularity build shapedb_container.sif docker-archive://shapedb_container.tar
