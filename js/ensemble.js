/*
 * The ensemble page.
 */

"use strict";

var ensemble;
var state;
var latestRevision;

function removeChildren(parent) {
  while (parent.firstChild) {
    parent.removeChild(parent.firstChild);
  }
}

function updatePage(json) {
  state = json;

  // Set the version options for the user to select.
  //
  var versions = [];
  for (let revision in state.versions) {
    versions.push(revision);
  }

  // wasteful comparison search but should not be an issue here
  versions.sort((a, b) => {
    const aa = Math.floor(a);
    const bb = Math.floor(b);
    if (aa < bb) { return -1; }
    if (aa > bb) { return 1; }
    return 0;
    }).reverse();

  var parent = document.getElementById('version');
  removeChildren(parent);
  for (let i = 0; i < versions.length; i++) {
    const el = document.createElement('option');
    el.value = versions[i];
    if (i === 0) { el.selected = true; }
    el.innerHTML = versions[i];

    parent.appendChild(el);
  }

  if (versions.length === 1) {
    parent.disabled = true;
  } else {
    parent.addEventListener("change",
			    (e) => { setVersion(e.target.value); });
  }

  // default to the latest version
  setVersion(state.latest_version);
}


// Send a save message if the contents have changed.
//
function saveUserContent() {
  const newval = document.getElementById('usercontent').value;

  const info = state.versions[latestRevision];
  if (newval === info.usernotes.user) {
    return;
  }

  let httpRequest = new XMLHttpRequest();
  if (!httpRequest) {
      alert("Unable to create a XMLHttpRequest!");
      return;
  }

  // Add the spinner to the whole page
  //
  let body = document.getElementsByTagName("body")[0];
  let spinner = new Spinner(spinopts);
  spinner.spin(body);

  // Update the local state once the save has been made
  httpRequest.addEventListener("load", function() {
    info.usernotes.user = newval;
  });
  httpRequest.addEventListener("error", function() {
      alert("Unable to save data!");
  });
  httpRequest.addEventListener("abort", function() {
      alert("Unable to save data!");
  });
  httpRequest.addEventListener("loadend", function() {
      spinner.stop();
  });

  const store = {name: ensemble,
		 revision: latestRevision,
		 usernotes: newval};
  httpRequest.open('POST', '/save/ensemble');
  httpRequest.setRequestHeader("Content-Type", "application/json;charset=UTF-8");
  httpRequest.send(JSON.stringify(store));
  
}


// Switch the display to the given version; minimal error checking
//
// If not the latest then page is "read-only".
//
function setVersion(revision) {

  const isLatest = state.latest_version === revision;

  let info = state.versions[revision];
  if (typeof info === 'undefined') {
    alert("Internal error: unknown revision: [" + revision + "]");
    return;
  }

  latestRevision = revision;

  document.getElementById("nstacks").innerHTML = info.nstacks.toString();
  document.getElementById("nmasters").innerHTML = info.nmasters.toString();

  /* handle the master list */
  let parent = document.getElementById('summarydata');
  let all_hulls_have_useraction = true;
  removeChildren(parent);
  for (let i = 0; i < info.nmasters; i++) {
    let m = info.masters[i];
    let el = document.createElement("tr");

    var td;
    td = document.createElement("td");
    td.innerHTML = m.masterid;
    el.appendChild(td);

    td = document.createElement("td");
    td.innerHTML = m.ncpts.toString();
    el.appendChild(td);

    td = document.createElement("td");
    var a = document.createElement("a");
    a.className = 'hullreview';
    a.href = "/" + info.name + "/" + revision + "/" + m.masterid;
    if (isLatest) { a.innerHTML = "Review"; } else { a.innerHTML = "View"; }
    td.appendChild(a);
    el.appendChild(td);

    td = document.createElement("td");
    if (m.useraction === '') {
      td.className = "undecided";
      td.innerHTML = "NO DECISION";
      all_hulls_have_useraction = false;
    } else {
      td.innerHTML = m.useraction;
    }
    el.appendChild(td);

    parent.appendChild(el);
  }

  // Update the "user content" section
  const usernotes = document.getElementById("usercontent");

  // We always use the user version if set, even if it is empty.
  //
  if (info.usernotes.user !== null) {
    usernotes.value = info.usernotes.user;
  } else {
    usernotes.value = info.usernotes.proposed;
  }

  const usersave = document.getElementById("saveusercontent");

  // Update the image
  parent = document.getElementById('overviewimage');
  parent.src = "/img/" + info.name + "/field." + info.name +
               ".v" + revision + ".png";

  // Is this an actionable page or not?
  // a) are we the latest version
  // b) do all the hulls have a user action?
  //
  const finish = document.getElementById('final');
  if (isLatest) {
    finish.style.display = 'block';
    document.getElementById('finish')
            .disabled = !all_hulls_have_useraction;

    usersave.addEventListener("click",
			      (e) => { saveUserContent(); });

  } else {
    finish.style.disply = 'none';

    usersave.disabled = true;
  }
}

const spinopts = {
    color: "#FF0000", opacity: 0.4,
    radius: 20, length: 20, width: 7
};


// Input is the ensemble value for the page.
//
function initialize(ens) {

  ensemble = ens;

  let httpRequest = new XMLHttpRequest();
  if (!httpRequest) {
      alert("Unable to create a XMLHttpRequest!");
      return;
  }

  // Add the spinner to the whole page
  //
  let body = document.getElementsByTagName("body")[0];
  let spinner = new Spinner(spinopts);
  spinner.spin(body);

  httpRequest.addEventListener("load", function() {
    updatePage(httpRequest.response);
  });
  httpRequest.addEventListener("error", function() {
      alert("Unable to load data!");
  });
  httpRequest.addEventListener("abort", function() {
      alert("Unable to load data!");
  });
  httpRequest.addEventListener("loadend", function() {
      spinner.stop();
  });

  // Do I need to add a cache-busting identifier?
  httpRequest.open('GET', '/api/' + ensemble + '?' +
                   (new Date()).getTime());
  httpRequest.responseType = 'json';
  httpRequest.send();
}
