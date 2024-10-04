Navigate to folder containing Dockerfile

Build docker:
sudo docker build .

Run docker with ports 8043 published:
sudo docker run --name=omics --rm --detach -p 8043:8043 6d4e8127ea77

## where "6d4e8127ea77" is the image number. Find with <sudo docker image ls>

Check if it's running:
sudo docker ps

Clear space:
sudo docker system prune -f -a

Stop docker:
sudo docker stop c0d70d53eeba

## where "c0d70d53eeba" is the container ID. Find with <sudo docker ps>
