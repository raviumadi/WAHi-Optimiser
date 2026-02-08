---
layout: default
title: Download
permalink: /download/
nav_order: 1
---

# Download WAH-*i*-Optimiser

Select your operating system to download the latest packaged release.

The installers include everything required to run the application.  
**MATLAB Runtime will be installed automatically if it is not already present.**

---

<div style="max-width:520px; margin:2rem auto; padding:1.5rem; border:1px solid #ddd; border-radius:12px; background:#fafafa; text-align:center;">

  <h3 style="margin-top:0;">Choose your platform</h3>

  <select id="osSelect"
          style="width:100%; padding:10px; font-size:1rem; border-radius:6px; border:1px solid #ccc; margin:1rem 0;">
    <option value="">— Select operating system —</option>
    <option value="mac">macOS (Apple Silicon)</option>
    <option value="win">Windows (64-bit)</option>
  </select>

  <a id="downloadBtn"
     href="#"
     style="display:inline-block; padding:12px 22px; margin-top:10px;
            font-size:1rem; font-weight:600; color:white;
            background:#2f80ed; border-radius:8px;
            text-decoration:none; pointer-events:none; opacity:0.5;">
    Download
  </a>

</div>

---

## Notes

- The packaged installers are generated using **MATLAB Application Compiler**
- All platforms run the **same optimiser core**
- Appearance may differ slightly due to system fonts and screen resolution. (Recommended: [install the Lato font](https://fonts.google.com/specimen/Lato))
- The application window size is fixed to ensure layout consistency

---

## License

The software is released under **GPLv3**.  
Academic and non-commercial research use is fully supported.

---

<script>
(function () {
  const select = document.getElementById('osSelect');
  const btn = document.getElementById('downloadBtn');

  const links = {
    mac: "{{ '/downloads/mac/wahi_gui_web.zip' | relative_url }}",
    win: "{{ '/downloads/winx64/wahi_gui_web.exe.zip' | relative_url }}"
  };

  select.addEventListener('change', function () {
    const os = select.value;

    if (links[os]) {
      btn.href = links[os];
      btn.style.pointerEvents = 'auto';
      btn.style.opacity = '1';
    } else {
      btn.href = '#';
      btn.style.pointerEvents = 'none';
      btn.style.opacity = '0.5';
    }
  });
})();
</script>