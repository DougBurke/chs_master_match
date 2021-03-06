/*
 * functions for the index page
 */

"use strict";

var usernotes = "";

// Since the data tables are handled very similarly, try and make the
// handling somewhat generic.
//
function tablefields() { return ['todo', 'review', 'completed', 'finished']; }
function mktablestruct() {
    return {todo: undefined, review: undefined,
            completed: undefined, finished: undefined};
}
var datatable = mktablestruct();


function plural(x) { if (x === 1) { return ""; } else { return "s"; } }

function nens(n) { return n.toString() + " ensemble" + plural(n); }

function add_ensemble_row(parent, ens) {

  // Add in 
  //   ensemble name
  //   number todo
  //   number masters
  //   number stacks
  //   revision

  const tr = document.createElement("tr");
  // tr.className = "ensemble";  // remove styling just for now

  let td = document.createElement("td");
  const a = document.createElement("a");
  a.href = ens.name;
  a.target = "_blank";
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
  ***/

    /*
  td = document.createElement("td");
    if ((ens.lastmodified.user === null) ||
        (ens.lastmodified.user === '')) {
    td.innerHTML = ens.lastmodified.proposed;
  } else {
    td.innerHTML = ens.lastmodified.user;
  }
  tr.appendChild(td);
    */
    
  td = document.createElement("td");
  td.innerHTML = ens.revision; // should already be a string
  tr.appendChild(td);

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

// Try and save information about the tables that we would want
// preserved on page reload, such as: size of table, current page
// number, row order.
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
  
  // Perhaps should have just saved these info structures
  const store = {length: mktablestruct(),
                 page: mktablestruct(),
                 order: mktablestruct()};

    for (var field of tablefields()) {
      const info = datatable[field].page.info();
      store.length[field] = info.length;
      store.page[field] = info.page,
      store.order[field] = datatable[field].order();
  }
         
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

  document.getElementById("ntodos").innerHTML = nens(json.todos.length);
  document.getElementById("nreviews").innerHTML = nens(json.reviews.length);
  document.getElementById("ncompleted").innerHTML = nens(json.completed.length);
  document.getElementById("nfinished").innerHTML = nens(json.finished.length);

  const errElem = document.getElementById("errors");
  const nerrors = json.errors.length;
  if (nerrors === 0) {
      errElem.innerHTML = "No errors";
  } else {
      errElem.innerHTML = nens(nerrors) + " - " + json.errors.join(", ");
      errElem.style.color = "red";
  }

  // TODO: do we have to clear out these divs?
  var parent = document.getElementById("todo-table-body");
  for (var ens of json.todos) {  
    add_ensemble_row(parent, ens);
  }

  parent = document.getElementById("review-table-body");
  for (var ens of json.reviews) {  
    add_ensemble_row(parent, ens);
  }

  parent = document.getElementById("completed-table-body");
  for (var ens of json.completed) {  
    add_ensemble_row(parent, ens);
  }

  parent = document.getElementById("finished-table-body");
  for (var ens of json.finished) {  
    add_ensemble_row(parent, ens);
  }

  // initialise the data tables
  datatable.todo = $('#todo-table').DataTable();
  datatable.review = $('#review-table').DataTable();
  datatable.completed = $('#completed-table').DataTable();
  datatable.finished = $('#finished-table').DataTable();

  if (json.datatable) {
    if (json.datatable.length) {
      for (var n in json.datatable.length) {
        datatable[n].page.len(json.datatable.length[n]);
      }
    }

    if (json.datatable.page) {
      for (var n in json.datatable.page) {
        datatable[n].page(json.datatable.page[n]);
      }
    }

    if (json.datatable.order) {
      for (var n in json.datatable.order) {
        datatable[n].order(json.datatable.order[n]);
      }
    }

    // not 100% convinced I have the draw argument correct here.
    for (var n of tablefields()) {
      datatable[n].draw(false);
    }
  }

  // Only set up the handler after setting the page length!
    for (var n of tablefields()) {
    for (var event of ['length', 'order', 'page']) {
      datatable[n].on(event, () => { saveDatatableSettings(); });
    }
  }
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
