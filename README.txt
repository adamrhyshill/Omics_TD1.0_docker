Build using docker only - no need for docker-compose.
Navigate to dir containing docker file.

Build docker:
sudo docker build .

Check if it's running and make a note of IMAGE_ID:
sudo docker ps

Find Image_ID:
docker image ls

Add flag to show even stopped processes:
sudo docker ps -a

Run docker:
sudo docker run --name=omics --rm --detach IMAGE_ID

Clean up all containers / images:
sudo docker system prune -f -a

See all images and containers and disk usage:
sudo docker system df

System defaults to running on port 8043 on local host (127.0.0.1)
