#This script packages the compiled app for distribution.
#This part of the code courtesy of Lasonic. Fixes the Qt library path Problem on Mac
BIN_FILE=Parametor.app/Contents/MacOS/Parametor
QWT_ROOT=/usr/local/qwt-6.1.0-svn/
for P in `otool -L $BIN_FILE | awk '{print $1}'` 
do 
    if [[ "$P" == *//* ]] 
    then 
        PSLASH=$(echo $P | sed 's,//,/,g')
        install_name_tool -change $P $PSLASH $BIN_FILE
    fi 
done 
 
QTDIR=/Users/smorrell/Qt5.1.0/5.1.0/clang_64
for F in `find $QTDIR/lib $QTDIR/plugins $QTDIR/qml  -perm 755 -type f` 
do 
    for P in `otool -L $F | awk '{print $1}'`
	do   
	    if [[ "$P" == *//* ]] 
	    then 
	        PSLASH=$(echo $P | sed 's,//,/,g')
	        install_name_tool -change $P $PSLASH $F
	    fi 
	 done
done

#Fix the QWT Dependence Path
QWT_LINK_PATH=$QWT_ROOT qwt.framework/Versions/6/qwt
install_name_tool -change qwt.framework/Versions/6/qwt $QWT_LINK_PATH $BIN_FILE

macdeployqt Parametor.app -dmg

#Copy QWT Framework to App MAnually
cp -r $QWT_ROOT"lib/qwt.framework" $PWD/Parametor.app/Contents/Frameworks