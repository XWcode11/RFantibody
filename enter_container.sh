
# first run the container
# docker run --name rfantibody_xw --gpus all -v .:/home --memory 20g -it rfantibody:xw

# start rfantibody:xw

docker start rfantibody_xw

#enter the container
docker exec -it rfantibody_xw /bin/bash