version="2.0"
year="2026"
rm -rf SemiSim.app
jpackage --type app-image \
	--app-version $version \
	--copyright "Brandon Li ($year)" \
	--name SemiSim \
	--icon ../images/icon.icns \
	--input ../target/ \
	--main-class electrodynamics.SemiSim \
	--main-jar SemiSim-$version.jar