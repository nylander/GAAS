install: update_profiles
	
update_profiles:
	sed -i.orig -e "s#VIRTUAL_ENV=\"\$$HOME/GAAS\"#VIRTUAL_ENV=\"${PWD}\"#" profiles/activate_nbis_env
	sed -i.orig -e "s#VIRTUAL_ENV=\"\$$HOME/GAAS\"#VIRTUAL_ENV=\"${PWD}\"#" profiles/activate_rackham_env

clean:
	mv profiles/activate_nbis_env{.orig,}
	mv profiles/activate_rackham_env{.orig,}
