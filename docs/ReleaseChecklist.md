# Creating a new release version

1. change version number in CommonCMake.txt
2. create binaries for all platforms with vagrant and fix build bugs
3. verify that all system tests pass
4. verify that crux-toolkit cucumber tests pass
5. add information to changelog on master
6. commit and push to master on github
7. create a release on github and upload binaries, use the tag naming convention rel-<major>-<minor>[-<patch>], e.g. rel-3-01 or rel-3-01-02
8. create a branch on github with the version number
