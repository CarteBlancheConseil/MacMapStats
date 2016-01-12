# MacMapStats

MacMapStats is an Xcode project containing statistical C utilities for :
- handling matrices,
- calculating simple datasets measures (min , max , standard deviation , ...), 
- computing analysis (factor analysis, classification, ...),
- read / write datas.

**Compilation :**
You can open MacMapStats.xcodeproj with Xcode then compile, or if you prefer to use the command line tool xcodebuild, open a terminal window and type a command like *xcodebuild install -configuration Release -Scheme MacMapStats ARCHS ="i386" -sdk macosx10.8*.

Note: MacMap projects are currently built with Xcode 6 against Mac OS X 10.8 SDK.

**Installation location :**
Default location for MacMap is "/Library/Framework", but like any Mac OS X framework, you can copy MacMapStats.framework to any standard location. (See https://developer.apple.com/library/mac/documentation/MacOSX/Conceptual/BPFrameworks/Tasks/InstallingFrameworks.html for more information).

