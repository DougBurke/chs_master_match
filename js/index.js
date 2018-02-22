function plural(x) { if (x === 1) { return ""; } else { return "s"; } }

function nens(n) { return n.toString() + " ensemble" + plural(n); }

function updatePage(json) {
  const ntodos = json.todos.length;
  const nreviews = json.reviews.length;
  const ncompleted = json.completed.length;

  document.getElementById("ntodos").innerHTML = nens(ntodos);
  document.getElementById("nreviews").innerHTML = nens(nreviews);
  document.getElementById("ncompleted").innerHTML = nens(ncompleted);

  // TODO: do we have to clear out these divs?
  var parent = document.getElementById("todo");
  for (let i = 0; i < ntodos; i++) {
    let ens = json.todos[i];

    // TODO: need to work out the number of remaining hulls
    //       or just do away with this concept (but would be nice
    //       to know which ones you've worked on)
    //
    let el = document.createElement("a");
    el.className = "ensemble";
    el.href = ens['name'];
    el.innerHTML = ens['name'] + "<br>" +
        ens['nmasters'].toString() +
        " remaining of " +
        ens['nmasters'].toString() + ", " +
        ens['nstacks'].toString() + " stack" + plural(ens['nstacks']) +
        "<br>N/A";

    parent.appendChild(el);
  }

  parent = document.getElementById("review");
  for (let i = 0; i < nreviews; i++) {
    let ens = json.reviews[i];

    // TODO: need to work out the number of remaining hulls
    let el = document.createElement("a");
    el.className = "ensemble";
    el.href = ens['name'];
    el.innerHTML = ens['name'] + "<br>" +
        ens['nmasters'].toString() + " hulls, " +
        ens['nstacks'].toString() + " stack" + plural(ens['nstacks']);

    parent.appendChild(el);
  }

  parent = document.getElementById("completed");
  for (let i = 0; i < ncompleted; i++) {
    let ens = json.completed[i];

    // TODO: need to work out the number of remaining hulls
    let el = document.createElement("a");
    el.className = "ensemble";
    el.href = ens['name'];
    el.innerHTML = ens['name'] + "<br>" +
        ens['nmasters'].toString() + " hulls, " +
        ens['nstacks'].toString() + " stack" + plural(ens['nstacks']);

    parent.appendChild(el);
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
