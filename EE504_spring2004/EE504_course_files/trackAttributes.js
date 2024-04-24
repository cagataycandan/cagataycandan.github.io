var urlList = new Array();

urlList["OCWHome"]="D|E|"; // Control for ocw home pages
urlList["Allcourseslist"]="D|E|"; // Control for ocw course list
urlList["Global"]="D|E|"; // Control for other global pages
//urlList["<dept content group Name>"]="E|E|"; 
//urlList["<course content group name>"]="E|E|E|E"; 

// Syntax:
// Global pages = <Netraker Flag> | <Sitewise Flag>
// Department pages = <Netraker Flag> | <Sitewise Flag>
// Course pages = <CHP Netraker Flag> | <CHP Sitewise Flag> | <Section // 				Netraker Flag> | <Section Sitewise Flag>
// An entry is made only if default behavior for that page needs to be // changed
// By Default for all pages, Netraker is disabled
// By Default Sitewise is enabled
// Each entry is assumed to be complete. 

function IsAnalyticsEnabled(contentGroup, contentSubGroup) {

var flagValues = null; // don't know the exact syntax. Please check.
var flags = null;

	if (contentGroup == "Global")	{
		if (contentSubGroup == "index") {
flagValues = urlList["OCWHome"];
		}
		else if (contentSubGroup == "all-courses")
		{
flagValues = urlList["Allcourseslist"];
		}
		else // it is not OCW Home page OR all courses list page
		{
flagValues = urlList["Global"];
		}

		if (!flagValues)
			return true; // enabled by default
		else {
			flags = flagValues.split("|");
			return (flags[1] == 'E');
		}
	} // End Global
	
	flagValues = urlList[contentGroup];
	if (!flagValues)
		return true; // enabled by default
	flags = flagValues.split("|");

	// Check if this is a department page or a Course Home page
	if (contentSubGroup == "" || contentSubGroup == "CourseHome")
	{
		return (flags[1] == 'E');
	} // End Department or CHP
	else // It is a course section.
	{
		return (flags[3] == 'E');
	}
}

function ComputeSurveyFlags (contentGroup, contentSubGroup) {

var flagValues = null; // don't know the exact syntax. Please check.
var flags = null;

	if (contentGroup == "Global")	{
		if (contentSubGroup == "index") {
flagValues = urlList["OCWHome"];
		}
		else if (contentSubGroup == "all-courses")
		{
flagValues = urlList["Allcourseslist"];
		}
		else // it is not OCW Home page OR all courses list page
		{
flagValues = urlList["Global"];
		}

		if (!flagValues)
			return false; // disabled by default
		else {
			flags = flagValues.split("|");
			return (flags[0] == 'E');
		}
	} // End Global
	
	flagValues = urlList[contentGroup];
	if (!flagValues)
		return false; // disabled by default
	flags = flagValues.split("|");

	// Check if this is a department page or a Course Home page
	if (contentSubGroup == "" || contentSubGroup == "CourseHome")
	{
		return (flags[0] == 'E');
	} // End Department OR CHP
	else // It is a course section.
	{
		return (flags[2] == 'E');
	}
}

function IsSurveyEnabled () {

var i;
var tmpContentGroup = "";
var tmpSubContentGroup = "";

//To check this function is supported by Client Browser
    if (! document.getElementsByTagName)
       {
         return false;
       }


var metaTags = document.getElementsByTagName("meta");
for (i = 0; i < metaTags.length; i++){
		if (metaTags[i].getAttribute("name") == "group"){
		tmpContentGroup = metaTags[i].getAttribute("content");
		}

		if (metaTags[i].getAttribute("name") == "subgroup"){
		tmpSubContentGroup = metaTags[i].getAttribute("content");
		}
}

return ComputeSurveyFlags (tmpContentGroup, tmpSubContentGroup);
}
