export async function loadHTMLIncludes(page) {
    const headerPlaceholder = document.getElementById('header-placeholder');
    const footerPlaceholder = document.getElementById('footer-placeholder');
    if (headerPlaceholder) {
        try {
            const response = await fetch('/static/html/_header.html');
            if (response.ok) {
                headerPlaceholder.innerHTML = await response.text();
                document.querySelectorAll('nav a').forEach(link => link.classList.remove('active'));
                const nav = document.getElementById(page);
                if (nav) {
                    nav.classList.add('active');
                }
            } else {
                headerPlaceholder.innerHTML = '<p style="color:red; text-align:center;">Error loading header.</p>';
            }
        } catch (error) {
            headerPlaceholder.innerHTML = '<p style="color:red; text-align:center;">Error loading header.</p>';
        }
    }

    if (footerPlaceholder) {
        try {
            const response = await fetch('/static/html/_footer.html');
            if (response.ok) {
                footerPlaceholder.innerHTML = await response.text();
            } else {
                footerPlaceholder.innerHTML = '<p style="color:red; text-align:center;">Error loading footer.</p>';
            }
        } catch (error) {
            footerPlaceholder.innerHTML = '<p style="color:red; text-align:center;">Error loading footer.</p>';
        }
    }
}
