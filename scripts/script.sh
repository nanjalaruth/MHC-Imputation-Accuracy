#Building a docker image
docker build -t impimage .

#Building a docker container
sudo docker run -i -t impimage /bin/bash 

#calling docker images
docker images

#calling docker containers
docker ps

#Built a docker image using quay.io, synched it with the github repository.
#Ran the nextflow pipeline on the server

#To be able to identify the container name:
sudo docker ps

# To be able to copy a file from a container use:
 docker cp romantic_leavitt:/usr/local/bin/SNP2HLA_package_v1.0.3/MakeReference/MakeReference.csh .
 
#I edited the shebanh in the two scripts and the path to the beagle.jar, beagle2linkage.jar and the linkage2beagle.jar
