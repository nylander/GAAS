install:
	mkdir bin
	find ${PWD}/annotation/ -perm +111 -type f -exec ln -s {} bin/ \;
	sed -i.orig -e "s#VIRTUAL_ENV=\"\$HOME/GAAS\"#VIRTUAL_ENV=\"${PWD}\"#" profiles/activate_nbis_env
	sed -i.orig -e "s#VIRTUAL_ENV=\"\$HOME/GAAS\"#VIRTUAL_ENV=\"${PWD}\"#" profiles/activate_rackham_env

clean:
	rm -rf bin
	mv profiles/activate_nbis_env{.orig,}
	mv profiles/activate_rackham_env{.orig,}
