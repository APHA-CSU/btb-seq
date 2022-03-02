set -e

################## DEPENDENCIES ######################

sudo apt-get install -y awscli 

#    ca-certificates \
#       gnupg \
#          lsb-release

# Add Docker’s official GPG key
#curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
#echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
#	    $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# Install Docker Engine
#sudo apt-get -y update
#printf "8\n" sudo | apt-get -y install docker-ce docker-ce-cli containerd.io
#sudo usermod -a -G docker "$USER" 

# Install Amazon ECS container agent
#sudo sh -c "echo 'net.ipv4.conf.all.route_localnet = 1' >> /etc/sysctl.conf"
#sudo sysctl -p /etc/sysctl.conf
#sudo apt-get install iptables-persistent
#sudo iptables -t nat -A PREROUTING -p tcp -d 169.254.170.2 --dport 80 -j DNAT --to-destination 127.0.0.1:51679
#sudo iptables -t nat -A OUTPUT -d 169.254.170.2 -p tcp -m tcp --dport 80 -j REDIRECT --to-ports 51679
#sudo iptables -A INPUT -i eth0 -p tcp --dport 51678 -j DROP
#sudo sh -c 'iptables-save > /etc/iptables/rules.v4'
#sudo mkdir -p /etc/ecs && sudo touch /etc/ecs/ecs.config
#sudo echo -e "ECS_DATADIR=/data \n
#	     ECS_ENABLE_TASK_IAM_ROLE=true \n
#	     ECS_ENABLE_TASK_IAM_ROLE_NETWORK_HOST=true \n
#	     ECS_LOGFILE=/log/ecs-agent.log \n
#	     ECS_AVAILABLE_LOGGING_DRIVERS=["json-file","awslogs"] \n
#	     ECS_LOGLEVEL=info \n
#	     ECS_CLUSTER=default" > /etc/ecs.config

#docker run --name ecs-agent \
#	--detach=true \
#	--restart=on-failure:10 \
#	--volume=/var/run:/var/run \
#	--volume=/var/log/ecs/:/log \
#	--volume=/var/lib/ecs/data:/data \
#	--volume=/etc/ecs:/etc/ecs \
#	--net=host \
#	--env-file=/etc/ecs/ecs.config \
#	amazon/amazon-ecs-agent:latest
