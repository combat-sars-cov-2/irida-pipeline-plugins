# COMBAT-SARS-COV-2 IRIDA Pipeline Plugins

[![Build Status](https://travis-ci.org/COMBAT-SARS-COV-2/irida-pipeline-plugins.svg?branch=master)](https://travis-ci.org/COMBAT-SARS-COV-2/irida-pipeline-plugins)

This project contains [IRIDA](http://irida.ca) pipeline plugins, custom developed for the [COMBAT-SARS-COV-2](http://combatsarscov2.org) project.

## The List of Plugins

1. artic-nanopore

## Installing Plugin to your local IRIDA

To install any of the plugins to your local IRIDA installation, please do the following:

Download a plugin package version of your choice from the [RELEASE](https://github.com/COMBAT-SARS-COV-2/irida-pipeline-plugins/releases) page.

```bash
cd /etc/irida/plugins
mv $plugin-package.jar to ./
```
Restart IRIDA services.

