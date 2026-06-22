version="2.0"
year="2026"
rm -rf output_linux
jpackage --type app-image \
	--app-version $version \
	--copyright "Brandon Li ($year)" \
	--name SemiSim \
	--icon ../images/icon.png \
	--input ../target/ \
	--dest output_linux \
	--main-class electrodynamics.SemiSim \
	--main-jar SemiSim-$version.jar

sh steam/ContentBuilder/scripts/build_linux.sh