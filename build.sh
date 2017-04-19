#!/bin/sh
cd $TRAVIS_BUILD_DIR/pkg
sbt ++$TRAVIS_SCALA_VERSION package
