
window.addEventListener("DOMContentLoaded", (ev) => {
    const alpha = 'abcdefghijklmnopqrstuvwxyz';
    let footnotes = document.getElementById("footnote-section");
    let idx = 0;

    for (let x of document.getElementsByClassName("footnote")) {
        let d = document.createElement("div");
        d.className = "endnote";
        d.id = `footnote-${idx}`

        let a = document.createElement("a");
        a.href = `#footnote-ref-${idx}`;
        a.text = `[${alpha[idx]}]`;
        a.style.paddingRight = '1rem';
        let c = document.createElement("span");
        c.innerHTML = x.innerHTML;
        d.appendChild(a);
        d.appendChild(c);
        console.log(c);
        footnotes.appendChild(d);
        idx += 1;
    }
    idx = 0;
    for (let x of document.getElementsByClassName("footnote-ref")) {
        // x.outerHTML = `<a href="#footnote-${idx}>${idx}</a>`;
        console.log(x);
        let ele = document.createElement("a");
        ele.href = `#footnote-${idx}`;

        x.id = `footnote-ref-${idx}`;
        ele.textContent = '[' + alpha[idx] + ']';
        x.appendChild(ele);
        idx += 1;
    }

});