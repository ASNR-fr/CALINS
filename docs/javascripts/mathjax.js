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
  }
};

document$.subscribe(() => {
  // Wait for MathJax to be fully initialized
  if (window.MathJax) {
    if (MathJax.startup && MathJax.startup.promise) {
      MathJax.startup.promise
        .then(() => MathJax.typesetPromise())
        .catch((err) => console.log('MathJax error:', err));
    } else if (MathJax.typesetPromise) {
      // Fallback for older versions
      MathJax.typesetPromise()
        .catch((err) => console.log('MathJax error:', err));
    }
  }
})