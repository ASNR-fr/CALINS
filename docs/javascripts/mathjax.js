window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"], ["$", "$"]],
    displayMath: [["\\[", "\\]"], ["$$", "$$"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  },
  startup: {
    ready() {
      MathJax.startup.defaultReady();
      document$.subscribe(({ body }) => {
        MathJax.typesetPromise([body])
          .catch(err => console.log('MathJax error:', err));
      });
    }
  }
};