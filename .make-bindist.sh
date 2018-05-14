#! /bin/sh
### Make a binary distribution for Mac OS X.
## To be run after "make install".
if [ "$TRAVIS_BRANCH" == master ]; then
set -e
echo "BINDISTNAME=@PACKAGE@-@VERSION@-@build@.tar.gz; BINDISTPATTERN=@PACKAGE@-*-@build@.tar.gz; prefix=@prefix@; srcdir=@srcdir@; abs_srcdir=@abs_srcdir@" > .bindist-vars.in
./config.status --quiet --file=.bindist-vars
. .bindist-vars
rm -f .bindist-vars.in .bindist-vars
(cd "$prefix" && tar cf - .) | gzip > $BINDISTNAME
echo Created $BINDISTNAME
## UPLOAD to github repository Normaliz-bindist
# This uses a secret deploy key.
# The key pair was created using:
#   ssh-keygen -t dsa -f .travis_ci_bindist_deploy_key
# The public key was uploaded to https://github.com/mkoeppe/Normaliz-bindist/settings/keys
# The private key was encrypted using
#   travis encrypt-file .travis_ci_bindist_deploy_key
# and appears as .travis_ci_bindist_deploy_key.enc in the repository.
# This also uploaded decryption keys in the form of environment variables on Travis CI.
if test x"$encrypted_571bff183f47_key" != x -a x"$encrypted_571bff183f47_iv" != x -a ! -r $srcdir/.travis_ci_bindist_deploy_key ; then
    openssl aes-256-cbc -K $encrypted_571bff183f47_key -iv $encrypted_571bff183f47_iv -in .travis_ci_bindist_deploy_key.enc -out $srcdir/.travis_ci_bindist_deploy_key -d
fi
if test -r $srcdir/.travis_ci_bindist_deploy_key; then 
    echo "Deployment key exists, attempting to upload"
    chmod 0600 $srcdir/.travis_ci_bindist_deploy_key
    export GIT_SSH_COMMAND="ssh -i $abs_srcdir/.travis_ci_bindist_deploy_key -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no"
    rm -Rf Normaliz-bindist
    git clone --depth 1 git@github.com:mkoeppe/Normaliz-bindist.git
    rm -f Normaliz-bindist/$BINDISTPATTERN
    cp $BINDISTNAME Normaliz-bindist/
    msg="Automatic upload from Travis CI, ${TRAVIS_REPO_SLUG} job=${TRAVIS_JOB_NUMBER} branch=${TRAVIS_BRANCH}"
    if [[ -n ${TRAVIS_TAG} ]]; then
        msg="${msg} tag=${TRAVIS_TAG}"
    fi
    msg="${msg} commit=${TRAVIS_COMMIT}"
    if [[ -n ${TRAVIS_PULL_REQUEST} ]]; then
        msg="${msg} pull_request=${TRAVIS_PULL_REQUEST}"
    fi
    (cd Normaliz-bindist && git --version && git add --all && git commit -m "${msg}" && git push)
fi
set +e
fi
