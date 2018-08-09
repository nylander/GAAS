# Profiles

These profiles are meant to make it easy to set up the environment needed for WebApollo without corrupting your own settings

## NBIS production server environment

To activate the NBIS production server environment use the profile `activate_nbis_env`.
```bash
# Change GAAS_HOME to where the GAAS repository is cloned
GAAS_HOME=$HOME/GAAS
source $GAAS_HOME/profiles/activate_nbis_env
```
Your prompt should now indicate the environment is active.

To restore your previous environment just type:
```bash
deactivate
```

## Rackham user environment

To activate the Rackham user environment use the profile `activate_rackham_env`.
```bash
# Change GAAS_HOME to where the GAAS repository is cloned
GAAS_HOME=$HOME/GAAS
source $GAAS_HOME/profiles/activate_rackham_env
```
Your prompt should now indicate the environment is active.

To restore your previous environment just type:
```bash
deactivate
```

