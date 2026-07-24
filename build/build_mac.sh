version="2.1"
year="2026"
rm -rf output_mac/SemiSim.app
jpackage --type app-image \
	--app-version $version \
	--copyright "Brandon Li ($year)" \
	--name SemiSim \
	--icon ../images/icon.icns \
	--input ../target/ \
	--dest output_mac \
	--main-class electrodynamics.SemiSim \
	--main-jar SemiSim-$version.jar \
	--java-options -XX:-TieredCompilation

sh steam/ContentBuilder/scripts/build_mac.sh