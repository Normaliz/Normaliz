for block in hollow_tri*.gz
do
    zcat $block | ../Dropbox/git_master/source/chunk/chunk
done
