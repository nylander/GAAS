check:
	@printf "Checking for tool dependencies\n"
	@printf "Ruby is installed ... "
	@if command -v ruby >/dev/null 2>&1 ; then printf "\033[0;32myes\033[0m : "; ruby -v; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Groovy is installed ... "
	@if command -v groovy >/dev/null 2>&1 ; then printf "\033[0;32myes\033[0m : "; groovy -v; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Perl is installed ... "
	@if command -v perl >/dev/null 2>&1 ; then printf "\033[0;32myes\033[0m : "; perl -v; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Rscript is installed ... "
	@if command -v Rscript >/dev/null 2>&1 ; then printf "\033[0;32myes\033[0m : "; Rscript --version; else printf "\033[0;31mno\033[0m\n"; fi

install: update_profiles
	
update_profiles:
	sed -i.orig -e "s#VIRTUAL_ENV=\"\$$HOME/GAAS\"#VIRTUAL_ENV=\"${PWD}\"#" profiles/activate_nbis_env
	sed -i.orig -e "s#VIRTUAL_ENV=\"\$$HOME/GAAS\"#VIRTUAL_ENV=\"${PWD}\"#" profiles/activate_rackham_env

clean:
	mv profiles/activate_nbis_env{.orig,}
	mv profiles/activate_rackham_env{.orig,}

get_groovy:
	curl -s "https://get.sdkman.io" | bash && source "${HOME}/.sdkman/bin/sdkman-init.sh" && sdk install groovy
