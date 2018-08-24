/*
 * functions for the index page
 */

"use strict";

var usernotes = "";

// datatable handling could be made more generic
//
var datatable_todo = undefined;
var datatable_review = undefined;
var datatable_completed = undefined;


function plural(x) { if (x === 1) { return ""; } else { return "s"; } }

function nens(n) { return n.toString() + " ensemble" + plural(n); }

function add_ensemble_row(parent, ens) {

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
  const info_todo = datatable_todo.page.info();
  const info_review = datatable_review.page.info();
  const info_completed = datatable_completed.page.info();

  const order_todo = datatable_todo.order();
  const order_review = datatable_review.order();
  const order_completed = datatable_completed.order();

  const store = {length: {todo: info_todo.length,
			  review: info_review.length,
			  completed: info_completed.length
			 },
		 page: {todo: info_todo.page,
			review: info_review.page,
			completed: info_completed.page
		       },
		 order: {todo: order_todo,
			 review: order_review,
			 completed: order_completed
			}
		};
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

  parent = document.getElementById("completed-table-body");
  for (let i = 0; i < ncompleted; i++) {
    let ens = json.completed[i];
    add_ensemble_row(parent, ens);
  }

  // initialise the data tables
  datatable_todo = $('#todo-table').DataTable();
  datatable_review = $('#review-table').DataTable();
  datatable_completed = $('#completed-table').DataTable();

  if (json.datatable) {
    if (json.datatable.length) {
      if (json.datatable.length.todo) {
        datatable_todo.page.len(json.datatable.length.todo);
      }
      if (json.datatable.length.review) {
        datatable_review.page.len(json.datatable.length.review);
      }
      if (json.datatable.length.completed) {
        datatable_completed.page.len(json.datatable.length.completed);
      }
    }

    if (json.datatable.page) {
      if (json.datatable.page.todo) {
        datatable_todo.page(json.datatable.page.todo);
      }
      if (json.datatable.page.review) {
        datatable_review.page(json.datatable.page.review);
      }
      if (json.datatable.page.completed) {
        datatable_completed.page(json.datatable.page.completed);
      }
    }

    if (json.datatable.order) {
      if (json.datatable.order.todo) {
        datatable_todo.order(json.datatable.order.todo);
      }
      if (json.datatable.order.review) {
        datatable_review.order(json.datatable.order.review);
      }
      if (json.datatable.order.completed) {
        datatable_completed.order(json.datatable.order.completed);
      }
    }

    // not 100% convinced I have the draw argument correct here.
    datatable_todo.draw(false);
    datatable_review.draw(false);
    datatable_completed.draw(false);
  }

  // Only set up the handler after setting the page length!
  for (var tbl of [datatable_todo, datatable_review, datatable_completed]) {
    for (var event of ['length', 'order', 'page']) {
      tbl.on(event, () => { saveDatatableSettings(); });
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
