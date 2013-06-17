umount -t vboxsf /share-dir;
rm -rf /share-dir;

mkdir /share-dir;
mount -t vboxsf  share-dir /share-dir;

