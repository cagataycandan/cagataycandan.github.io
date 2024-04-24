// OCW custom code -- BEGIN
var i;
var MainContentGroup = "";
var SubContentGroup = "";

if (document.getElementsByTagName) 
{
var metaTags = document.getElementsByTagName("meta");
	
// Retrieve content group and content sub-groups from the metatags.
for (i = 0; i < metaTags.length; i++) {
if (metaTags[i].getAttribute("name") == "group"){
		MainContentGroup = metaTags[i].getAttribute("content");
	}

	if (metaTags[i].getAttribute("name") == "subgroup"){
		SubContentGroup = metaTags[i].getAttribute("content");
	}
}
}
// OCW Custom code -- END

var wtl_TagVer = 6;
var wtl_FWD = 0;
var wtl_url = document.URL;
var wtl_title = document.title;
var wtl_TagID = 144556;
var wtl_SID = "144556";
var wtl_Offset = "-500";
WTL_TAG = new Image;
WTL_TAG.ID = "WTL_TAG";

var ORDER= "";
var SERVER= "";
var INVOICE= "";
var CARTVIEW= "";
var CARTADD= "";
var CARTREMOVE= "";
var CHECKOUT= "";
var CARTBUY= "";
var ADCAMPAIGN= "";

//if (IsAnalyticsEnabled(MainContentGroup, SubContentGroup)) {
//alert ("sitewise is enabled");
var contentGroupAndSubgroup = MainContentGroup + "," + SubContentGroup;
wtl_Tag6(wtl_TagID,wtl_SID,wtl_url,wtl_title,contentGroupAndSubgroup);
//}

// NOTE:
// The above is the "final" state implementation. For the time being
// until we get comfortable with Sitewise, we should remove the if
// condition. The code will look like the following:

//var contentGroupAndSubgroup = MainContentGroup + "," + SubContentGroup;
//wtl_Tag6(wtl_TagID,wtl_SID,wtl_url,wtl_title,contentGroupAndSubgroup);
