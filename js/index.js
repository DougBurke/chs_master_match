/*
 * functions for the index page
 */

"use strict";

var usernotes = "";

var datatable_todo = undefined;
var datatable_review = undefined;


function plural(x) { if (x === 1) { return ""; } else { return "s"; } }

function nens(n) { return n.toString() + " ensemble" + plural(n); }

function add_ensemble_row(parent, ens) {

  const tr = document.createElement("tr");
  // tr.className = "ensemble";  // remove styling just for now

  let td = document.createElement("td");
  const a = document.createElement("a");
  a.href = ens.name;
  a.innerHTML = ens.name;
  td.appendChild(a);
  tr.appendChild(td);

  // We consider those hulls in the masters array which have
  // useraction.user set to a non-empty value to be "done".
  //
  let done = 0;
  for (var h of ens.masters) {
      const ua = h.useraction.user;
      if ((ua !== null) && (ua !== "")) { done += 1; }
  }

  const rem = ens.nmasters - done;

  td = document.createElement("td");
  td.innerHTML = rem.toString();
  tr.appendChild(td);

  td = document.createElement("td");
  td.innerHTML = ens.nmasters.toString();
  tr.appendChild(td);

  td = document.createElement("td");
  td.innerHTML = ens.nstacks.toString();
  tr.appendChild(td);

  /***
      At the moment the last-modified field only refers to
      any ensemble-level changes (i.e. user comment)
      and does not include changes made to each hull
      so it is probably confusing

  td = document.createElement("td");
  if (ens.lastmodified.user === null) {
    td.innerHTML = ens.lastmodified.proposed;
  } else {
    td.innerHTML = ens.lastmodified.user;
  }
  tr.appendChild(td);
  ***/

  parent.appendChild(tr);
}

// Send a save message if the contents have changed.
//
function saveUserContent() {
  const newval = document.getElementById('usercontent').value;

  if (newval === usernotes) {
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
    usernotes = newval;
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

  const store = {usernotes: newval};
  httpRequest.open('POST', '/save/summary');
  httpRequest.setRequestHeader("Content-Type", "application/json;charset=UTF-8");
  httpRequest.send(JSON.stringify(store));
  
}

// For now this means that the display length of one of the tables
// has been changed. Could also store current page number.
//
function saveDatatableSettings() {

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

  const store = {length: {todo: datatable_todo.page.len(),
			  review: datatable_review.page.len()}};
  httpRequest.open('POST', '/save/datatable');
  httpRequest.setRequestHeader("Content-Type", "application/json;charset=UTF-8");
  httpRequest.send(JSON.stringify(store));
  
}

function updatePage(json) {

  document.getElementById("usercontent").value = json.usernotes;
  usernotes = json.usernotes;

  // Allow the user to add notes
  document.getElementById("saveusercontent")
      .addEventListener("click",
			(e) => { saveUserContent(); });

  const ntodos = json.todos.length;
  const nreviews = json.reviews.length;
  const ncompleted = json.completed.length;

  document.getElementById("ntodos").innerHTML = nens(ntodos);
  document.getElementById("nreviews").innerHTML = nens(nreviews);
  document.getElementById("ncompleted").innerHTML = nens(ncompleted);

  // TODO: do we have to clear out these divs?
  var parent = document.getElementById("todo-table-body");
  for (let i = 0; i < ntodos; i++) {
    let ens = json.todos[i];

    // TODO: need to work out the number of remaining hulls
    //       or just do away with this concept (but would be nice
    //       to know which ones you've worked on)
    //
    add_ensemble_row(parent, ens);
  }

  parent = document.getElementById("review-table-body");
  for (let i = 0; i < nreviews; i++) {
    let ens = json.reviews[i];
    add_ensemble_row(parent, ens);
  }

  /*** not used at present
  parent = document.getElementById("completed-table-body");
  for (let i = 0; i < ncompleted; i++) {
    let ens = json.completed[i];
    add_ensemble_row(parent, ens);
  }
  ***/

  // initialise the data tables
  datatable_todo = $('#todo-table').DataTable();
  datatable_review = $('#review-table').DataTable();

  if (json.datatable && json.datatable.length) {
      if (json.datatable.length.todo) {
	datatable_todo.page.len(json.datatable.length.todo).draw();
      }
      if (json.datatable.length.review) {
	datatable_review.page.len(json.datatable.length.review).draw();
      }
  }

  // Only set up the handler after setting the page length!
  datatable_todo.on('length', () => {
    saveDatatableSettings();
  }); 
  datatable_review.on('length', () => {
    saveDatatableSettings();
  }); 
}

const spinopts = {
    color: "#FF0000", opacity: 0.4,
    radius: 20, length: 20, width: 7
};

function initialize() {
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
  httpRequest.open('GET', '/api/summary?' +
                   (new Date()).getTime());
  httpRequest.responseType = 'json';
  httpRequest.send();
}
