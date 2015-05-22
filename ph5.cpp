#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "FireLog.h"
#include "FireUtils.hpp"
#include "version.h"
#include "ph5.h"

using namespace std;
using namespace ph5;

static void help() {
    cout << "ph5 Pythagorean Hodograph Quintic curve library" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << endl;
    cout << "Copyright 2015, Karl Lew" << endl;
    cout << "https://github.com/ph5/README.md" << endl;
}

bool parseArgs(int argc, char *argv[]) {
    firelog_level(FIRELOG_INFO);

    if (argc <= 1) {
        return false;
    }

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == 0) {
            // empty argument
        } else if (strcmp("-version", argv[i]) == 0) {
            cout << "{\"version\":\"" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << "\"}" << endl;
            return false;
        } else if (strcmp("-warn", argv[i]) == 0) {
            firelog_level(FIRELOG_WARN);
        } else if (strcmp("-error", argv[i]) == 0) {
            firelog_level(FIRELOG_ERROR);
        } else if (strcmp("-info", argv[i]) == 0) {
            firelog_level(FIRELOG_INFO);
        } else if (strcmp("-debug", argv[i]) == 0) {
            firelog_level(FIRELOG_DEBUG);
        } else if (strcmp("-trace", argv[i]) == 0) {
            firelog_level(FIRELOG_TRACE);
        } else {
            LOGERROR1("unknown ph5 argument: '%s'", argv[i]);
            return false;
        }
    }
    return true;
}

int main(int argc, char *argv[]) {
	cout << "PH5 hello" << endl;
    bool argsOk = parseArgs(argc, argv);
    if (!argsOk) {
        help();
		return -1;
    }

	/*
    Mat image;
    if (imagePath) {
        LOGTRACE1("Reading image: %s", imagePath);
        image = imread(imagePath);
        if (!image.data) {
            LOGERROR1("main() imread(%s) failed", imagePath);
            exit(-1);
        }
    } else {
        LOGDEBUG("No image specified.");
    }

    switch (uimode) {
    case UI_STILL:
        uiStill(pipelinePath.c_str(), image, argMap, isTime, jsonIndent);
        break;
    case UI_VIDEO:
        uiVideo(pipelinePath.c_str(), argMap);
        break;
    default:
        LOGERROR("Unknown UI mode");
        exit(-1);
    }

    if (outputPath) {
        if (!imwrite(outputPath, image)) {
            LOGERROR1("Could not write image to: %s", outputPath);
            exit(-1);
        }
        LOGTRACE1("Image written to: %s", outputPath);
    }

*/
    return 0;
}
