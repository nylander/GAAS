check:
	@printf "Checking for tool dependencies\n"
	@printf "Ruby is installed ... "
	@if command -v ruby >/dev/null 2>&1 ; then printf "\033[0;32myes\033[0m : "; ruby -v; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Groovy is installed ... "
	@if command -v groovy >/dev/null 2>&1 ; then printf "\033[0;32myes\033[0m : "; groovy -v; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Perl is installed ... "
	@if command -v perl >/dev/null 2>&1 ; then printf "\033[0;32myes\033[0m : "; perl -v; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Perl lib Clone is installed ... "
	@if perl -MClone -e 1  2>/dev/null; then printf "\033[0;32myes\033[0m \n"; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Perl lib Moose is installed ... "
	@if perl -MMoose -e 1 2>/dev/null; then printf "\033[0;32myes\033[0m \n"; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Perl lib Graph::Directed is installed ..."
	@if perl -MGraph::Directed -e 1 2>/dev/null; then printf "\033[0;32myes\033[0m \n"; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Perl lib LWP::UserAgent is installed ..."
	@if perl -MLWP::UserAgent -e 1 2>/dev/null; then printf "\033[0;32myes\033[0m \n"; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Perl lib Statistics::R is installed ... " 
	@if perl -MStatistics::R -e 1 2>/dev/null; then printf "\033[0;32myes\033[0m \n"; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Perl lib JSON is installed ... " 
	@if perl -MJSON -e 1 2>/dev/null; then printf "\033[0;32myes\033[0m \n"; else printf "\033[0;31mno\033[0m\n"; fi
	@printf "Rscript is installed ... "
	@if command -v Rscript >/dev/null 2>&1 ; then printf "\033[0;32myes\033[0m : "; Rscript --version; else printf "\033[0;31mno\033[0m\n"; fi

install: update_profiles
	
update_profiles:
	sed -i.orig -e "s#VIRTUAL_ENV=\"\$$HOME/GAAS\"#VIRTUAL_ENV=\"${PWD}\"#" profiles/activate_nbis_env
	sed -i.orig -e "s#VIRTUAL_ENV=\"\$$HOME/GAAS\"#VIRTUAL_ENV=\"${PWD}\"#" profiles/activate_rackham_env
	sed -i.orig -e "s#VIRTUAL_ENV=\"\$$HOME/GAAS\"#VIRTUAL_ENV=\"${PWD}\"#" profiles/activate_env
clean:
	mv profiles/activate_nbis_env{.orig,}
	mv profiles/activate_rackham_env{.orig,}
	mv profiles/activate_local_env{.orig,}

get_groovy:
	curl -s "https://get.sdkman.io" | bash && source "${HOME}/.sdkman/bin/sdkman-init.sh" && sdk install groovy
