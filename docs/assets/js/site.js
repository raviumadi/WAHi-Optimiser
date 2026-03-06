document.addEventListener("DOMContentLoaded", () => {
  const root = document.documentElement;
  const storedTheme = localStorage.getItem("wahi-theme");
  const preferredDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
  const initialTheme = storedTheme || (preferredDark ? "dark" : "light");

  function applyTheme(theme) {
    root.setAttribute("data-theme", theme);
    const label = document.querySelector("[data-theme-label]");
    if (label) {
      label.textContent = theme === "dark" ? "Switch to light theme" : "Switch to dark theme";
    }
  }

  applyTheme(initialTheme);

  const themeToggle = document.getElementById("themeToggle");
  if (themeToggle) {
    themeToggle.addEventListener("click", () => {
      const nextTheme = root.getAttribute("data-theme") === "dark" ? "light" : "dark";
      localStorage.setItem("wahi-theme", nextTheme);
      applyTheme(nextTheme);
    });
  }

  const navToggle = document.getElementById("siteNavToggle");
  const nav = document.getElementById("siteNav");
  if (navToggle && nav) {
    navToggle.addEventListener("click", () => {
      const isOpen = nav.classList.toggle("is-open");
      navToggle.setAttribute("aria-expanded", String(isOpen));
    });

    window.addEventListener("resize", () => {
      if (window.innerWidth > 820) {
        nav.classList.remove("is-open");
        navToggle.setAttribute("aria-expanded", "false");
      }
    });
  }
});
